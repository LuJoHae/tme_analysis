"""Dataset providers for single-cell, bulk iAtlas, direct papers, spatial, and TCGA."""

from .bulk_iatlas import IATLAS_COHORTS, load_iatlas_cohort
from .bulk_papers import PAPER_DATASETS, load_genentech_egad, load_paper_h5ad
from .single_cell import (
    load_jerby_arnon,
    load_ma_liver,
    load_maynard,
    load_sade_feldman,
)
from .spatial import compute_spatial_graph, load_spatial_dataset
from .tcga import load_tcga_project

__all__ = [
    "load_sade_feldman",
    "load_jerby_arnon",
    "load_maynard",
    "load_ma_liver",
    "IATLAS_COHORTS",
    "load_iatlas_cohort",
    "PAPER_DATASETS",
    "load_paper_h5ad",
    "load_genentech_egad",
    "load_spatial_dataset",
    "compute_spatial_graph",
    "load_tcga_project",
]
