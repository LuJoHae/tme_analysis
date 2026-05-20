import pandas as pd
import numpy as np
import pyensembl
import re
import mygene


def _normalize_gene(id, ensembl, original_id=None):
    if id is None:
        return pd.Series({"original.id": original_id})
    try:
        gene = ensembl.gene_by_id(id.split(".")[0])
        gene_normalized = pd.Series(gene.to_dict())
        if (gene_normalized["gene_name"] is None) or (gene_normalized["gene_name"] == ""):
            gene_normalized["gene_name"] = "UNKNOWN"
        gene_normalized["EnsemblRelease.release"] = gene_normalized["genome"].release
        gene_normalized["EnsemblRelease.species"] = gene_normalized["genome"].species.latin_name
        gene_normalized = gene_normalized.drop("genome")
        gene_normalized["original.id"] = original_id
        return gene_normalized
    except:
        return pd.Series({"original.id": original_id})


def _detect_pre_id_transform(var_names, ensembl):
    if np.all([bool(re.fullmatch(r"ENSG\d{11}", name)) for name in var_names]):
        pre_id_transform = None
    elif np.all([bool(re.fullmatch(r"ENST\d{11}", name)) for name in var_names]):
        pre_id_transform = "transcript"
    elif np.all([bool(re.fullmatch(r"\d{1,10}", name)) for name in var_names]):
        pre_id_transform = "entrez"
    elif sum(name in ensembl.gene_names() for name in var_names) >= 0.1 * len(var_names):
        pre_id_transform = "hugo"
    else:
        raise ValueError("Could not automatically determine pre_id_transform!")
    return pre_id_transform


def _normalize_all_genes(var_names, pre_id_transform=None):
    ensembl = pyensembl.EnsemblRelease(release=111, species="human")
    original_ids = list(var_names)

    if pre_id_transform == "auto":
        pre_id_transform = _detect_pre_id_transform(var_names, ensembl)

    match pre_id_transform:
        case "transcript":
            ids = list()
            for id in list(var_names):
                try:
                    gene_id = ensembl.transcript_by_id(id).gene_id
                    ids.append(gene_id)
                except ValueError:
                    ids.append(None)
        case "entrez":
            result = mygene.MyGeneInfo().getgenes(list(var_names), fields="symbol,name,ensembl")
            ids = list()
            for gene in result:
                if "ensembl" not in gene.keys():
                    ids.append(None)
                    continue
                if type(gene["ensembl"]) is list:
                    ids.append(None)
                    continue
                ids.append(gene["ensembl"]["gene"])
        case "hugo":
            ids = []
            for gene_name in var_names:
                try:
                    ensembl_ids = ensembl.gene_ids_of_gene_name(gene_name)
                except ValueError:
                    ids.append(None)
                    continue
                if len(ensembl_ids) == 1:
                    ids.append(ensembl_ids[0])
                else:
                    ids.append(None)
        case None:
            ids = list(var_names)
    genes = pd.concat([_normalize_gene(gene_id, ensembl, original_id) for gene_id, original_id in zip(ids, original_ids)], axis=1).T
    return genes


def _set_dtypes(var_matrix):
    var_matrix["contig"] = var_matrix["contig"].astype("category")
    var_matrix["start"] = var_matrix["start"].astype("string")
    var_matrix["end"] = var_matrix["end"].astype("string")
    var_matrix["strand"] = var_matrix["strand"].astype("category")
    var_matrix["biotype"] = var_matrix["biotype"].astype("category")
    var_matrix["gene_id"] = var_matrix["gene_id"].astype("string")
    var_matrix["gene_name"] = var_matrix["gene_name"].astype("string")
    var_matrix["EnsemblRelease.release"] = var_matrix["EnsemblRelease.release"].astype("category")
    var_matrix["EnsemblRelease.species"] = var_matrix["EnsemblRelease.species"].astype("category")
    var_matrix["original.id"] = var_matrix["original.id"].astype("string")
    return var_matrix


def norm_genes(adata, pre_id_transform):
    var_matrix = _normalize_all_genes(adata.var_names, pre_id_transform=pre_id_transform)
    assert np.all(var_matrix["original.id"] == adata.var.index)
    var_matrix = var_matrix.dropna()
    var_matrix = _set_dtypes(var_matrix)
    adata = adata[:, var_matrix.index]
    adata.var = var_matrix
    adata.var.set_index("gene_id", inplace=True)
    adata = adata[:, sorted(list(adata.var.index))]
    return adata
