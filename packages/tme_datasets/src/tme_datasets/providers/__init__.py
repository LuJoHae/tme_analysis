"""Dataset providers for single-cell, bulk iAtlas, direct papers, spatial, and TCGA."""

from .bulk_iatlas import IATLAS_COHORTS, load_iatlas_cohort
from .bulk_papers import PAPER_DATASETS, load_genentech_egad, load_paper_h5ad
from .single_cell import (
    load_gse179994,
    load_jerby_arnon,
    load_ma_liver,
    load_maynard,
    load_sade_feldman,
    load_yost,
)
from .atlas_single_cell import (
    load_azizi_brca,
    load_becker_coad,
    load_biermann_brainmet,
    load_borcherding_ccrcc,
    load_cheng_pancancer,
    load_durante_uvm,
    load_khaliq_cc,
    load_kim_luad,
    load_leader_nsclc,
    load_lu_hcc,
    load_pelka_crc,
    load_pu_ptc,
    load_qian_pancancer,
    load_sharma_hcc,
    load_vazquez_ov,
    load_zhang_myeloid,
    load_zhang_tnbc,
)
from .spatial import compute_spatial_graph, load_spatial_dataset
from .tcga import load_tcga_project

__all__ = [
    "load_sade_feldman",
    "load_jerby_arnon",
    "load_maynard",
    "load_ma_liver",
    "load_yost",
    "load_gse179994",
    "load_pelka_crc",
    "load_azizi_brca",
    "load_qian_pancancer",
    "load_cheng_pancancer",
    "load_leader_nsclc",
    "load_kim_luad",
    "load_becker_coad",
    "load_khaliq_cc",
    "load_borcherding_ccrcc",
    "load_sharma_hcc",
    "load_lu_hcc",
    "load_pu_ptc",
    "load_durante_uvm",
    "load_biermann_brainmet",
    "load_vazquez_ov",
    "load_zhang_tnbc",
    "load_zhang_myeloid",
    "IATLAS_COHORTS",
    "load_iatlas_cohort",
    "PAPER_DATASETS",
    "load_paper_h5ad",
    "load_genentech_egad",
    "load_spatial_dataset",
    "compute_spatial_graph",
    "load_tcga_project",
]

