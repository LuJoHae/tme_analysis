import logging
import os
import shutil
import tarfile
import urllib
from pathlib import Path
import re
from typing import Any, Mapping, Sequence
import anndata as ad
import mygene
import numpy as np
import pandas as pd
import polars as pl
import pyensembl
import scipy.sparse as sp

logger = logging.getLogger("gene_utils")


def ensure_ensembl_release_installed(
    release: int = 111,
    species: str = "human",
    ensembl_dir: Path | None = None,
) -> pyensembl.EnsemblRelease:
    """Ensure the specified Ensembl release is downloaded and indexed by pyensembl.

    Delegates to tme_datasets.preprocessing.gene_normalization.
    """
    from tme_datasets.preprocessing.gene_normalization import (
        ensure_ensembl_release_installed as _impl,
    )

    return _impl(release=release, species=species, ensembl_dir=ensembl_dir)


def normalize_genes_to_ensembl(
    adata: ad.AnnData,
    release: int = 111,
    species: str = "human",
    ensembl_dir: Path | None = None,
    drop_unmapped: bool = True,
    aggregation: str = "sum",
) -> ad.AnnData:
    """Convert AnnData var_names to canonical Ensembl gene IDs and enrich .var attributes.

    Delegates to tme_datasets.preprocessing.gene_normalization.
    """
    from tme_datasets.preprocessing.gene_normalization import (
        normalize_genes_to_ensembl as _impl,
    )

    return _impl(
        adata=adata,
        release=release,
        species=species,
        ensembl_dir=ensembl_dir,
        drop_unmapped=drop_unmapped,
        aggregation=aggregation,
    )


def norm_genes(
    adata: ad.AnnData,
    pre_id_transform: str | None = None,
    release: int = 111,
) -> ad.AnnData:
    """Normalize AnnData gene identifiers to Ensembl release (backward-compatible)."""
    return normalize_genes_to_ensembl(adata, release=release, drop_unmapped=True)


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


def rpkm_to_tpm(rpkm_df: pd.DataFrame) -> pd.DataFrame:
    """Normalize each column such that its sum is strictly 10^6."""
    return rpkm_df.div(rpkm_df.sum(axis=0), axis=1) * 1e6