import shutil
import tarfile
import urllib
from pathlib import Path
import pandas as pd
import numpy as np
import anndata as ad
import pyensembl
import re
import mygene


def rpkm_to_tpm(rpkm_array: np.ndarray) -> np.ndarray:
    """
    Converts an array of RPKM values to TPM.
    """
    # Sum of all RPKM values in the sample
    rpkm_sum = np.sum(rpkm_array, axis=0)

    # Scale each RPKM value by the sum, then multiply by 10^6
    tpm_array = (rpkm_array / rpkm_sum) * 1e6

    return tpm_array


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

    # Handle duplicate gene_ids that can arise when multiple input IDs
    # (e.g., Entrez IDs) map to the same Ensembl gene. We aggregate
    # expression values across duplicates by summing, and keep the first
    # annotation row for each gene_id.
    if not adata.var.index.is_unique:
        gene_ids = adata.var.index.to_numpy()
        unique_ids, inverse = np.unique(gene_ids, return_inverse=True)

        # Aggregate X by summing columns belonging to the same gene_id
        X = adata.X
        # Build a sparse-like aggregation using pandas for robustness
        X_dense = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        X_df = pd.DataFrame(X_dense, index=adata.obs_names, columns=gene_ids)
        X_agg = X_df.T.groupby(level=0).sum().T  # sum duplicates
        # Reorder columns to match unique_ids order
        X_agg = X_agg.loc[:, unique_ids]

        # Keep first occurrence of var annotations per unique gene_id
        var_dedup = adata.var[~adata.var.index.duplicated(keep="first")]
        var_dedup = var_dedup.loc[unique_ids]

        adata = ad.AnnData(
            X=X_agg.to_numpy(),
            obs=adata.obs.copy(),
            var=var_dedup.copy(),
        )

    adata = adata[:, sorted(list(adata.var.index))]
    return adata


def read_gmt(file_path):
    """
    Parses a GMT file into a dictionary of lists.

    Returns:
        dict: Keys are gene set names, values are lists of gene symbols.
    """
    gene_sets = {}

    with open(file_path, 'r') as f:
        for line in f:
            # Strip whitespace and split by tab
            parts = line.strip().split('\t')

            # parts[0] is the name, parts[1] is the description
            name = parts[0]
            genes = parts[2:]

            gene_sets[name] = genes

    return gene_sets


def gmt_to_long_df(file_path):
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            name = parts[0]
            for gene in parts[2:]:
                data.append({'gene_set': name, 'gene': gene})

    return pd.DataFrame(data)


def download(url: str, dest: Path) -> None:
    """Download `url` to `dest` using a browser-like User-Agent header.

    Required because some servers (incl. science.bostongene.com) reject
    requests that come with the default Python-urllib User-Agent (HTTP 403).
    """
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": (
                "Mozilla/5.0 (X11; Linux x86_64) "
                "AppleWebKit/537.36 (KHTML, like Gecko) "
                "Chrome/120.0.0.0 Safari/537.36"
            ),
            "Accept": "*/*",
        },
    )
    with urllib.request.urlopen(request) as response, open(dest, "wb") as out_file:
        shutil.copyfileobj(response, out_file)


class GeneSet(object):
    def __init__(self, name, descr, genes):
        self.name = name
        self.descr = descr
        self.genes = set(genes)
        self.genes_ordered = list(genes)

    def __str__(self):
        return '{}\t{}\t{}'.format(self.name, self.descr, ', '.join(self.genes))

    def __repr__(self):
        return '{}\t{}\t{}'.format(self.name, self.descr, ', '.join(self.genes))


def read_gene_sets(gmt_file):
    """
    Return dict {geneset_name : GeneSet object}

    :param gmt_file: str, path to .gmt file
    :return: dict
    """
    gene_sets = {}
    with open(gmt_file) as handle:
        for line in handle:
            items = line.strip().split('\t')
            name = items[0].strip()
            description = items[1].strip()
            genes = set([gene.strip() for gene in items[2:]])
            gene_sets[name] = GeneSet(name, description, genes)

    return gene_sets


def ssgsea_score(ranks, genes):
    """According to bagaev"""
    common_genes = list(set(genes).intersection(set(ranks.index)))
    if not len(common_genes):
        return pd.Series([0] * len(ranks.columns), index=ranks.columns)
    sranks = ranks.loc[common_genes]
    return (sranks ** 1.25).sum() / (sranks ** 0.25).sum() - (
                len(ranks.index) - len(common_genes) + 1) / 2


def ssgsea_formula(data, gene_sets, rank_method='max'):
    """
    Return DataFrame with ssgsea scores
    Only overlapping genes will be analyzed

    :param data: pd.DataFrame, DataFrame with samples in columns and variables in rows
    :param gene_sets: dict, keys - processes, values - bioreactor.gsea.GeneSet
    :param rank_method: str, 'min' or 'max'.
    :return: pd.DataFrame, ssgsea scores, index - genesets, columns - patients
    """

    ranks = data.T.rank(method=rank_method, na_option='bottom')

    return pd.DataFrame({gs_name: ssgsea_score(ranks, gene_sets[gs_name].genes)
                         for gs_name in list(gene_sets.keys())})


def median_scale(data, clip=None):
    mad = (data - data.mean()).abs().mean()
    c_data = (data - data.median()) / mad
    if clip is not None:
        return c_data.clip(-clip, clip)
    return c_data


def download_github_file(owner, repo, commit_hash, file_path, output_filename):
    url = f"https://raw.githubusercontent.com/{owner}/{repo}/{commit_hash}/{file_path}"

    try:
        urllib.request.urlretrieve(url, output_filename)
        print(f"Success: File saved to {output_filename}")
    except urllib.error.HTTPError as e:
        print(f"HTTP Error: {e.code} - {e.reason}. Check if the file exists at this commit.")
    except Exception as e:
        print(f"An error occurred: {e}")


def download_from_cbioportal(filename: str, output_dir: Path) -> None:
    tar_filename = filename
    archive_path = output_dir / tar_filename
    extract_path = output_dir / tar_filename.removesuffix(".tar.gz")
    download(f"https://datahub.assets.cbioportal.org/{filename}", archive_path)
    with tarfile.open(archive_path, "r:gz") as tar:
        tar.extractall(path=extract_path)


def print_tree(node: Path, prefix: str = "") -> None:
    """Performs a depth-first traversal to print the tree topology."""
    if not node.is_dir():
        return

    children = list(node.iterdir())
    # Sort children to group directories and files consistently
    children.sort(key=lambda x: (x.is_file(), x.name))

    for index, child in enumerate(children):
        is_last = index == len(children) - 1
        connector = "└── " if is_last else "├── "

        print(f"{prefix}{connector}{child.name}")

        if child.is_dir():
            extension = "    " if is_last else "│   "
            print_tree(child, prefix + extension)


def calculate_maf_tmb(df: pd.DataFrame, capture_size_mb: float = 30.0, vaf_threshold: float = 0.05) -> pd.DataFrame:
    """
    Computes Tumor Mutational Burden (TMB) across a cohort from a MAF-formatted DataFrame.
    """
    # 1. Compute VAF (handling potential division by zero)
    df['VAF'] = np.where(df['t_depth'] > 0, df['t_alt_count'] / df['t_depth'], 0)

    # 2. Define the qualifying domain S (protein-altering variants in coding regions)
    qualifying_classifications = {
        'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation',
        'Frame_Shift_Ins', 'Frame_Shift_Del', 'In_Frame_Ins',
        'In_Frame_Del', 'Translation_Start_Site', 'Splice_Site'
    }

    # 3. Isolate valid somatic mutations
    qualifying_subset = df[
        (df['Variant_Classification'].isin(qualifying_classifications)) &
        (df['VAF'] >= vaf_threshold)
    ]

    # 4. Aggregate |S| per sample
    mutation_counts = qualifying_subset.groupby('Tumor_Sample_Barcode').size().reset_index(name='Qualifying_Mutations')

    # 5. Compute the scalar projection TMB = |S| / K
    mutation_counts['TMB_Score'] = mutation_counts['Qualifying_Mutations'] / capture_size_mb

    # 6. Reintegrate samples with |S| = 0 (those dropped during the filtering step)
    all_samples = pd.DataFrame({'Tumor_Sample_Barcode': df['Tumor_Sample_Barcode'].unique()})
    tmb_matrix = pd.merge(all_samples, mutation_counts, on='Tumor_Sample_Barcode', how='left')

    # Fill NaN values with 0 for samples with no qualifying mutations
    tmb_matrix['Qualifying_Mutations'] = tmb_matrix['Qualifying_Mutations'].fillna(0).astype(int)
    tmb_matrix['TMB_Score'] = tmb_matrix['TMB_Score'].fillna(0.0)

    return tmb_matrix