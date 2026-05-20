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


class DatasetSingleCell(datalair.Dataset):
    """Datalair Dataset class for all single cell datasets."""

    def __init__(self) -> None:
        """Initialize this dataset class as a datalair.Dataset class with namespace `DatasetSingleCell`."""
        super().__init__(namespace="DatasetSingleCell")


def _read_azizi(ds, lair):
    lair.safe_derive(ds, overwrite=False)
    adatas_matrix = list()
    adatas_counts = list()
    for filepath in sorted(list(lair.get_dataset_filepaths(ds).values())):
        if filepath.name[0] == ".":
            continue
        adata = ad.read_h5ad(filepath)
        adata.X = csr_matrix(adata.X)
        if "matrix" in filepath.name:
            adatas_matrix.append(adata)
        elif "counts" in filepath.name:
            adatas_counts.append(adata)
        else:
            raise ValueError(f"Unknown file type: {filepath}")
    adata_matrix = ad.concat(adatas_matrix, axis=0, join="outer")
    adata_counts = ad.concat(adatas_counts, axis=0, join="outer")

    adata_matrix = norm_genes(adata_matrix, pre_id_transform=None)
    adata_counts = norm_genes(adata_counts, pre_id_transform="hugo")

    adata_counts.obs["original.barcode"] = adata_counts.obs.index
    adata_matrix.var.drop(columns=["original.id"], inplace=True)
    adata_counts.var.drop(columns=["original.id"], inplace=True)

    adata = ad.concat([adata_matrix, adata_counts], join="outer", merge="first")
    adata = adata[:, sorted(list(adata.var.index))]
    return adata


def _read_pelka(ds, lair):
    lair.safe_derive(ds, overwrite=False)
    adata = ad.read_h5ad(lair.get_dataset_filepaths(ds)["adata.h5ad"])
    print(adata.shape)
    sc.pp.filter_genes(adata, min_counts=1)
    adata.var["gene_ids"] = [gene_id.split(".")[0] for gene_id in adata.var["gene_ids"] if gene_id[:4]=="ENSG"]
    adata.var.set_index("gene_ids", inplace=True)
    adata.obs["original.barcode"] = adata.obs.index
    adata.obs.reset_index(inplace=True, drop=True)
    adata = norm_genes(adata, pre_id_transform="auto")
    print(adata.shape)
    return adata


def _standard_read(ds, lair):
    lair.safe_derive(ds, overwrite=False)
    adata = ad.read_h5ad(lair.get_dataset_filepaths(ds)["adata.h5ad"])
    print(adata.shape)
    adata.obs["original.barcode"] = adata.obs.index
    adata.obs.reset_index(inplace=True, drop=True)
    adata = norm_genes(adata, pre_id_transform="auto")
    print(adata.shape)
    return adata


class SingleCellDataProcessStep01(DatasetSingleCell):
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
        ad.settings.allow_write_nullable_strings = True
        datasets = [
            singlecellrnasignature.adata.AziziSingleCellMapDiverse2018Adata(),
            singlecellrnasignature.adata.BeckerSinglecellAnalysesDefine2022Adata(),
            # singlecellrnasignature.adata.BiTumorImmuneReprogramming2021Adata(), # Derivation is broken! Raw data not downloaded? Broken!
            singlecellrnasignature.adata.BiermannDissectingTreatmentnaiveEcosystem2022Adata(),
            singlecellrnasignature.adata.BorcherdingMappingImmuneEnvironment2021Adata(),
            singlecellrnasignature.adata.ChengPancancerSinglecellTranscriptional2021Adata(),
            singlecellrnasignature.adata.DuranteSinglecellAnalysisReveals2020Adata(),
            singlecellrnasignature.adata.JerbyArnonCancerCellProgram2018Adata(),
            singlecellrnasignature.adata.KhaliqRefiningColorectalCancer2022Adata(),
            singlecellrnasignature.adata.KimSinglecellRNASequencing2020Adata(),
            # singlecellrnasignature.adata.KrishnaSinglecellSequencingLinks2021Adata(), # Need to manually convert Seurat to h5ad file
            singlecellrnasignature.adata.LeaderSinglecellAnalysisHuman2021Adata(),
            singlecellrnasignature.adata.LuSinglecellAtlasMulticellular2022Adata(),
            singlecellrnasignature.adata.PelkaSpatiallyOrganizedMulticellular2021Adata(),
            singlecellrnasignature.adata.PuSinglecellTranscriptomicAnalysis2021Adata(),
            singlecellrnasignature.adata.QianPancancerBlueprintHeterogeneous2020aAdata(),
            singlecellrnasignature.adata.SharmaOncofetalReprogrammingEndothelial2020Adata(),
            # Index of adata doubled
            singlecellrnasignature.adata.VazquezOvarianCancerMutational2022Adata(),
            singlecellrnasignature.adata.ZhangSinglecellAnalysesReveal2021Adata(),
            singlecellrnasignature.adata.ZhangSinglecellAnalysisReveals2022Adata()
            # yes, this is a different study; just happens that there are two people called Zhang that work on the same research and both don't write creative titles :)
        ]
        adatas = dict()
        for ds in datasets:
            print(ds._dataset_name)
            lair.safe_derive(ds)
            match ds._dataset_name:
                case "AziziSingleCellMapDiverse2018Adata":
                    adata = _read_azizi(ds, lair)
                case "PelkaSpatiallyOrganizedMulticellular2021Adata":
                    adata = _read_pelka(ds, lair)
                case _:
                    adata = _standard_read(ds, lair)
            adata = adata[adata.X.sum(
                axis=1) > 100, :]  # only load what could be a cell and not unfiltered junk
            print(type(adata.X))
            adata.write_h5ad(
                output_dir.joinpath(f"{ds._dataset_name}.h5ad"))
            adatas[ds._dataset_name] = adata


class SingleCellDataProcessStep02(DatasetSingleCell):
    storage_path = Path("/storage/halu").resolve()

    def derive(self, lair: datalair.Lair) -> None:
        output_dir = lair.get_path(self)
        adatas = dict()
        filelist = sorted(list(
            lair.get_dataset_filepaths(SingleCellDataProcessStep01()).values()
        ))
        for filepath in filelist:
            print(filepath.stem)
            adata = ad.read_h5ad(filepath, backed="r+")
            adata.obs.reset_index(inplace=True, drop=True)
            # adata = adata[np.random.choice(adata.obs_names, 100, replace=True), :]
            adata = adata.to_memory()
            adata.obs.reset_index(inplace=True, drop=True)
            adatas[filepath.stem] = adata

        cancer_types = {
            singlecellrnasignature.adata.AziziSingleCellMapDiverse2018Adata()._dataset_name:
                (["geo_id", "patient", "tissue", "original.barcode"],
                 ["geo_id", "patient", "organ", "original.barcode"], "BRCA",
                 "primary tumor"),
            singlecellrnasignature.adata.BeckerSinglecellAnalysesDefine2022Adata()._dataset_name:
                (["sample", "geo_id", "original.barcode"],
                 ["sample", "geo_id", "original.barcode"], "COAD", "primary tumor"),
            singlecellrnasignature.adata.BiermannDissectingTreatmentnaiveEcosystem2022Adata()._dataset_name:
                (["patient", "batch", "organ", "cell_type_main", "original.barcode"],
                 ["patient", "batch", "organ", "cell_type", "original.barcode"], "SKCM",
                 "brain metastasis"),
            singlecellrnasignature.adata.BorcherdingMappingImmuneEnvironment2021Adata()._dataset_name:
                (["original.barcode"], ["original.barcode"], "ccRCC", "primary tumor"),
            singlecellrnasignature.adata.ChengPancancerSinglecellTranscriptional2021Adata()._dataset_name:
                (["patient", "batch", "tissue", "cancer", "cancer_type",
                  "original.barcode"],
                 ["patient", "batch", "organ", "cancer", "cancer_type",
                  "original.barcode"], "Pan", "primary tumor"),
            # only myeloid cells excell sheets in supp data
            singlecellrnasignature.adata.DuranteSinglecellAnalysisReveals2020Adata()._dataset_name:
                (["sample", "geo_id", "original.barcode"],
                 ["sample", "geo_id", "original.barcode"], "UVM", "primary tumor"),
            singlecellrnasignature.adata.JerbyArnonCancerCellProgram2018Adata()._dataset_name:
                (["samples", "Cohort", "cell.types"],
                 ["sample", "cohort", "cell_types"], "SKCM", "primary tumor"),
            singlecellrnasignature.adata.KhaliqRefiningColorectalCancer2022Adata()._dataset_name:
                (["samples", "Condition"], ["sample", "is_tumor"], "CC",
                 "primary tumor"),
            singlecellrnasignature.adata.KimSinglecellRNASequencing2020Adata()._dataset_name:
                (["original.barcode"], ["original.barcode"], "LUAD", "primary tumor"),
            singlecellrnasignature.adata.LeaderSinglecellAnalysisHuman2021Adata()._dataset_name:
                (["0", "batch"], ["original.barcode", "batch"], "NSCLC",
                 "primary tumor"),
            singlecellrnasignature.adata.LuSinglecellAtlasMulticellular2022Adata()._dataset_name:
                (["sample", "patient", "site", "celltype", "original.barcode"],
                 ["sample", "patient", "organ", "cell_type", "original.barcode"], "HCC",
                 "primary tumor and metastasis"),
            singlecellrnasignature.adata.PelkaSpatiallyOrganizedMulticellular2021Adata()._dataset_name:
                (["SPECIMEN_TYPE", "PROCESSING_TYPE", "PatientTypeID",
                  "HistologicTypeSimple"],
                 ["is_tumor", "cell_type", "patient", "histology"], "CRC",
                 "primary tumor and metastasis"),
            singlecellrnasignature.adata.PuSinglecellTranscriptomicAnalysis2021Adata()._dataset_name:
                (["patient", "biopsy_site", "geo_id", "original.barcode"],
                 ["patient", "organ", "geo_id", "original.barcode"], "PTC",
                 "primary tumor"),  # papillary thyroid carcinoma
            singlecellrnasignature.adata.QianPancancerBlueprintHeterogeneous2020aAdata()._dataset_name:
                (["Source Name", "Characteristics[individual]",
                  "Characteristics[organism part]", "Characteristics[disease]",
                  "original.barcode"],
                 ["sample", "patient", "organ", "cancer_type", "original.barcode"],
                 "Pan", "primary_tumor"),  # lung, colorectal, ovary and breast cancer
            singlecellrnasignature.adata.SharmaOncofetalReprogrammingEndothelial2020Adata()._dataset_name:
                (["barcode"], ["original.barcode"], "HCC", "primary tumor"),
            singlecellrnasignature.adata.VazquezOvarianCancerMutational2022Adata()._dataset_name:
                (["sample", "patient_id", "cell_type", "tumor_site",
                  "original.barcode"],
                 ["sample", "patient", "cell_type", "organ", "original.barcode"], "OV",
                 "primary tumor"),
            singlecellrnasignature.adata.ZhangSinglecellAnalysesReveal2021Adata()._dataset_name:
                (["original.barcode"], ["original.barcode"], "TNBR", "primary tumor"),
            # triple negative
            singlecellrnasignature.adata.ZhangSinglecellAnalysisReveals2022Adata()._dataset_name:
                (["geo_id", "patient_id"], ["geo_id", "patient"], "SKCM",
                 "primary tumor")
            # not totally correct, there is also another skin cancer type in the dataset
        }

        print(20 * "=")
        for name, adata in adatas.items():
            print(name)
            adata.obs = adata.obs[cancer_types[name][0]]
            adata.obs.columns = cancer_types[name][1]
            adata.obs["cancer_code"] = cancer_types[name][2]
            adata.obs["staging"] = cancer_types[name][3]
            adata.obs["dataset"] = name
            adata.var["ensembl_id"] = adata.var.index
            if "original.id" in adata.var.columns:
                adata.var = adata.var.drop(columns="original.id")
            adata.write_h5ad(output_dir.joinpath(f"{name}.h5ad"))


class SingleCellDataProcessStep03(DatasetSingleCell):
    storage_path = Path("/storage/halu").resolve()

    def derive(self, lair: datalair.Lair) -> None:
        ad.settings.allow_write_nullable_strings = True
        output_dir = lair.get_path(self)
        filelist = sorted(list(
            lair.get_dataset_filepaths(SingleCellDataProcessStep02()).values()
        ))
        for filepath in filelist:
            name = filepath.stem
            adata = ad.read_h5ad(filepath, backed="r+")
            adata.obs.reset_index(inplace=True, drop=True)
            # adata = adata[np.random.choice(adata.obs_names, 100, replace=True), :]
            adata = adata.to_memory()
            adata.obs.reset_index(inplace=True, drop=True)
            sc.pp.filter_cells(adata, min_genes=300)
            adata.write_h5ad(output_dir.joinpath(f"{name}.h5ad"))


class SingleCellDataProcessStep04(DatasetSingleCell):
    storage_path = Path("/storage/halu").resolve()

    def derive(self, lair: datalair.Lair) -> None:
        ad.settings.allow_write_nullable_strings = True
        output_dir = lair.get_path(self)
        adatas = list()
        names = list()
        filelist = sorted(list(
            lair.get_dataset_filepaths(SingleCellDataProcessStep03()).values()
        ))
        for filepath in filelist:
            print(filepath.stem)
            adata = ad.read_h5ad(filepath, backed="r")
            adata.obs["dataset"] = filepath.stem
            adatas.append(adata)
            names.append(filepath.stem)

        gene_ids = sorted(list(set.union(*[set(adata.var_names) for adata in adatas])))
        xs = [csr_matrix((adata.n_obs, len(gene_ids)), dtype=np.int64) for adata in
              adatas]
        bdatas = [ad.AnnData(X=x, var=pd.DataFrame(index=gene_ids), obs=adata.obs) for
                  (x, adata) in zip(xs, adatas)]
        for i, (adata, bdata, name) in enumerate(zip(adatas, bdatas, names)):
            print(i)
            cadata = bdata.copy()
            adata_mem = adata.to_memory()
            adata_mem.X = csr_matrix(adata_mem.X, dtype=np.int64)

            # Map adata_mem's var_names to column positions in bdata
            cdata_var_index = {name: idx for idx, name in enumerate(cadata.var_names)}
            col_map = np.array([cdata_var_index[v] for v in adata_mem.var_names],
                               dtype=np.int64)

            # Remap columns from adata_mem's CSR into bdata's column space
            src = adata_mem.X.tocsr()
            new_indices = col_map[src.indices]  # remap column indices
            cadata.X = csr_matrix(
                (src.data, new_indices, src.indptr),
                shape=cadata.shape,
                dtype=np.int64,
            )
            cadata.write_h5ad(output_dir.joinpath(f"{name}.h5ad"))
            del adata_mem


class SingleCellDataProcessStep05(DatasetSingleCell):
    storage_path = Path("/storage/halu").resolve()

    def derive(self, lair: datalair.Lair) -> None:
        ad.settings.allow_write_nullable_strings = True
        output_dir = lair.get_path(self)
        ds = single_cell_datasets.SingleCellDataProcessStep04()
        filter_opts = {name: {"min_genes": min_genes, "max_genes": max_genes,
                              "min_counts": min_counts, "max_counts": max_counts,
                              "mt_content": mt_content} for
                       name, min_genes, max_genes, min_counts, max_counts, mt_content in
                       zip(
                           ["PuSinglecellTranscriptomicAnalysis2021Adata.h5ad",
                            "BiermannDissectingTreatmentnaiveEcosystem2022Adata.h5ad",
                            "KimSinglecellRNASequencing2020Adata.h5ad",
                            "ChengPancancerSinglecellTranscriptional2021Adata.h5ad",
                            "KhaliqRefiningColorectalCancer2022Adata.h5ad",
                            "LuSinglecellAtlasMulticellular2022Adata.h5ad",
                            "LeaderSinglecellAnalysisHuman2021Adata.h5ad",
                            "DuranteSinglecellAnalysisReveals2020Adata.h5ad",
                            "AziziSingleCellMapDiverse2018Adata.h5ad",
                            "QianPancancerBlueprintHeterogeneous2020aAdata.h5ad",
                            "BorcherdingMappingImmuneEnvironment2021Adata.h5ad",
                            "SharmaOncofetalReprogrammingEndothelial2020Adata.h5ad",
                            "PelkaSpatiallyOrganizedMulticellular2021Adata.h5ad",
                            "JerbyArnonCancerCellProgram2018Adata.h5ad",
                            "BeckerSinglecellAnalysesDefine2022Adata.h5ad"],
                           [350, 400, 400, 400, 300, 300, 300, 300, 300, 300, 400, 300,
                            300, 300, 300],
                           [3_000, 4_500, 2_500, 1_400, 3_000, 3_000, 3_000, 3_500,
                            3_000, 3_500, 3_000, 2_500, 3_000, 4_500, 2_500],
                           [500, 300, 300, 500, 300, 300, 300, 300, 300, 300, 300, 300,
                            300, 300, 300],
                           [10_000, 15_000, 15_000, 2_000, 14_000, 15_000, 14_000,
                            15_000, 9_000, 15_000, 15_000, 10_000, 15_000, 20_000,
                            5_000],
                           [10, 10, 12, 5, 20, 15, 25, 30, 20, 20, 13, 15, 20, 20, 5]
                           , strict=True)}

        for name, filepath in lair.get_dataset_filepaths(ds).items():
            print(name)
            adata = ad.read_h5ad(filepath)
            bdata = adata[:1, :].copy()
            bdata = gene_utils.norm_genes(bdata, pre_id_transform=None)
            adata.var = bdata.var
            adata.var["mt"] = (adata.var["contig"] == "MT").values
            sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                                       log1p=False, inplace=True)
            sc.pp.filter_cells(adata, min_genes=filter_opts[name]["min_genes"])
            sc.pp.filter_cells(adata, min_counts=filter_opts[name]["min_counts"])
            sc.pp.filter_cells(adata, max_genes=filter_opts[name]["max_genes"])
            sc.pp.filter_cells(adata, max_counts=filter_opts[name]["max_counts"])
            adata = adata[
                adata.obs.pct_counts_mt < filter_opts[name]["mt_content"], :].copy()
            adata.write_h5ad(output_dir.joinpath(name))


class SingleCellDataProcessStep06(DatasetSingleCell):
    storage_path = Path("/storage/halu").resolve()

    def derive(self, lair: datalair.Lair) -> None:
        ad.settings.allow_write_nullable_strings = True
        output_dir = lair.get_path(self)
        ds = single_cell_datasets.SingleCellDataProcessStep05()
        infiles = list(lair.get_dataset_filepaths(ds).values())
        ad.experimental.concat_on_disk(
            infiles, output_dir.joinpath("adata.h5ad"), join="inner", axis=0
        )


class SingleCellDataProcessStep07(DatasetSingleCell):
    storage_path = Path("/storage/halu").resolve()

    def derive(self, lair: datalair.Lair) -> None:
        ad.settings.allow_write_nullable_strings = True
        output_dir = lair.get_path(self)
        ds = SingleCellDataProcessStep06()
        lair.get_dataset_filepaths(ds)
        adata = ad.read_h5ad(lair.get_dataset_filepaths(ds)["adata.h5ad"])
        adata = gene_utils.norm_genes(adata, pre_id_transform=None)
        adata.obs.reset_index(inplace=True)
        junk_genes_masks = (adata.var["contig"] == "MT") | (
                    adata.var["contig"].apply(len) >= 3)
        adata = adata[:, ~junk_genes_masks]
        adata.write_h5ad(output_dir.joinpath("adata.h5ad"))
