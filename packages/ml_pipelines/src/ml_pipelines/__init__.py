"""Gene utils to be used with anndata."""

from __future__ import annotations

__version__ = "0.1.0"

from gene_utils._gene_utils import (
    norm_genes,
    download,
    ssgsea_score,
    ssgsea_formula,
    median_scale,
    read_gene_sets,
    read_gmt,
    download_github_file,
    download_from_cbioportal,
)

__all__ = [
    "__version__",
    "norm_genes",
    "download",
    "ssgsea_score",
    "ssgsea_formula",
    "median_scale",
    "read_gene_sets",
    "read_gmt",
    "download_github_file",
    "download_from_cbioportal"
]