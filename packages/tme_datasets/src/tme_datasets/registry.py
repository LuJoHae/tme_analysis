"""Declarative registry of all single-cell, spatial, bulk, and signature datasets."""

from __future__ import annotations

from typing import Mapping
from returns.maybe import Maybe, Nothing, Some

from .models import DatasetSpec
from .types import Modality

DATASET_REGISTRY: Mapping[str, DatasetSpec] = {
    # -------------------------------------------------------------------------
    # Single-Cell Reference Pool (docs/dataset_catalog/01_single_cell_reference_datasets.md)
    # -------------------------------------------------------------------------
    "GSE120575": DatasetSpec(
        id="GSE120575",
        title="Sade-Feldman et al. 2018 Melanoma T-cell ICB Response",
        modality=Modality.SINGLE_CELL,
        cancer_type="Melanoma",
        platform="Smart-seq2",
        organ=Some("Skin"),
        has_response_labels=True,
        n_samples_or_cells=Some(16288),
    ),
    "GSE115978": DatasetSpec(
        id="GSE115978",
        title="Jerby-Arnon et al. 2018 Melanoma T-cell Exclusion",
        modality=Modality.SINGLE_CELL,
        cancer_type="Melanoma",
        platform="Smart-seq2",
        organ=Some("Skin"),
        has_response_labels=False,
        n_samples_or_cells=Some(7186),
    ),
    "GSE125449": DatasetSpec(
        id="GSE125449",
        title="Ma et al. 2019 Liver Cancer (HCC / ICC) TME",
        modality=Modality.SINGLE_CELL,
        cancer_type="Liver Cancer",
        platform="10x_v2",
        organ=Some("Liver"),
        has_response_labels=False,
        n_samples_or_cells=Some(5115),
    ),
    "GSE123813": DatasetSpec(
        id="GSE123813",
        title="Yost et al. 2019 BCC Clonal Replacement",
        modality=Modality.SINGLE_CELL,
        cancer_type="Basal Cell Carcinoma",
        platform="10x_5prime",
        organ=Some("Skin"),
        has_response_labels=False,
        n_samples_or_cells=Some(3500),
    ),
    "Maynard_NSCLC": DatasetSpec(
        id="Maynard_NSCLC",
        title="Maynard et al. 2020 Therapy-Induced NSCLC Evolution",
        modality=Modality.SINGLE_CELL,
        cancer_type="Non-Small Cell Lung",
        platform="10x_v3",
        organ=Some("Lung"),
        has_response_labels=False,
        n_samples_or_cells=Some(3000),
    ),
    "GSE179994": DatasetSpec(
        id="GSE179994",
        title="Tietscher et al. Pan-Cancer T-cell Infiltrate Atlas",
        modality=Modality.SINGLE_CELL,
        cancer_type="Pan-Cancer",
        platform="10x_v3",
        organ=Some("Multi-organ"),
        has_response_labels=False,
        n_samples_or_cells=Some(150849),
    ),

    # -------------------------------------------------------------------------
    # Bulk iAtlas Validation Cohorts (docs/dataset_catalog/02_bulk_validation_cohorts.md)
    # -------------------------------------------------------------------------
    "Hugo-iAtlas": DatasetSpec(
        id="Hugo-iAtlas",
        title="Hugo et al. 2016 Melanoma anti-PD-1",
        modality=Modality.BULK_RNA,
        cancer_type="Melanoma",
        platform="RNA-seq",
        organ=Some("Skin"),
        has_response_labels=True,
        n_samples_or_cells=Some(27),
    ),
    "Riaz-iAtlas": DatasetSpec(
        id="Riaz-iAtlas",
        title="Riaz et al. 2017 Melanoma Nivolumab",
        modality=Modality.BULK_RNA,
        cancer_type="Melanoma",
        platform="RNA-seq",
        organ=Some("Skin"),
        has_response_labels=True,
        n_samples_or_cells=Some(107),
    ),
    "Liu-iAtlas": DatasetSpec(
        id="Liu-iAtlas",
        title="Liu et al. 2019 Melanoma anti-PD-1",
        modality=Modality.BULK_RNA,
        cancer_type="Melanoma",
        platform="RNA-seq",
        organ=Some("Skin"),
        has_response_labels=True,
        n_samples_or_cells=Some(122),
    ),
    "Gide-iAtlas": DatasetSpec(
        id="Gide-iAtlas",
        title="Gide et al. 2019 Melanoma anti-PD-1 +/- anti-CTLA-4",
        modality=Modality.BULK_RNA,
        cancer_type="Melanoma",
        platform="RNA-seq",
        organ=Some("Skin"),
        has_response_labels=True,
        n_samples_or_cells=Some(91),
    ),
    "Rosenberg-iAtlas": DatasetSpec(
        id="Rosenberg-iAtlas",
        title="Rosenberg et al. 2016 / IMvigor210 Bladder Atezolizumab",
        modality=Modality.BULK_RNA,
        cancer_type="Urothelial Bladder",
        platform="RNA-seq",
        organ=Some("Bladder"),
        has_response_labels=True,
        n_samples_or_cells=Some(347),
    ),
    "Padron-iAtlas": DatasetSpec(
        id="Padron-iAtlas",
        title="Padron et al. 2022 Pancreatic anti-PD-1 + CD40",
        modality=Modality.BULK_RNA,
        cancer_type="Pancreatic Ductal",
        platform="RNA-seq",
        organ=Some("Pancreas"),
        has_response_labels=True,
        n_samples_or_cells=Some(93),
    ),
    "Anders-iAtlas": DatasetSpec(
        id="Anders-iAtlas",
        title="Anders et al. 2021 Triple-Negative Breast Atezolizumab",
        modality=Modality.BULK_RNA,
        cancer_type="Triple-Negative Breast",
        platform="RNA-seq",
        organ=Some("Breast"),
        has_response_labels=True,
        n_samples_or_cells=Some(31),
    ),
    "McDermott-iAtlas": DatasetSpec(
        id="McDermott-iAtlas",
        title="McDermott et al. 2018 Renal Cell Atezolizumab +/- Bevacizumab",
        modality=Modality.BULK_RNA,
        cancer_type="Clear Cell Renal Cell",
        platform="RNA-seq",
        organ=Some("Kidney"),
        has_response_labels=True,
        n_samples_or_cells=Some(263),
    ),
    "Choueiri-iAtlas": DatasetSpec(
        id="Choueiri-iAtlas",
        title="Choueiri et al. 2016 Renal Cell Nivolumab",
        modality=Modality.BULK_RNA,
        cancer_type="Clear Cell Renal Cell",
        platform="RNA-seq",
        organ=Some("Kidney"),
        has_response_labels=True,
        n_samples_or_cells=Some(16),
    ),

    # -------------------------------------------------------------------------
    # Direct Paper ICI Response Cohorts (dataset_papers/README.md)
    # -------------------------------------------------------------------------
    "Auslander": DatasetSpec(
        id="Auslander",
        title="Auslander et al. 2018 Melanoma Response",
        modality=Modality.BULK_RNA,
        cancer_type="Melanoma",
        platform="RNA-seq",
        organ=Some("Skin"),
        has_response_labels=True,
    ),
    "Chen-CTLA4": DatasetSpec(
        id="Chen-CTLA4",
        title="Chen et al. 2016 Melanoma anti-CTLA-4",
        modality=Modality.BULK_RNA,
        cancer_type="Melanoma",
        platform="RNA-seq",
        organ=Some("Skin"),
        has_response_labels=True,
    ),
    "Chen-PD1": DatasetSpec(
        id="Chen-PD1",
        title="Chen et al. 2016 Melanoma anti-PD-1",
        modality=Modality.BULK_RNA,
        cancer_type="Melanoma",
        platform="RNA-seq",
        organ=Some("Skin"),
        has_response_labels=True,
    ),
    "Freeman": DatasetSpec(
        id="Freeman",
        title="Freeman et al. 2022 Melanoma ICB Response",
        modality=Modality.BULK_RNA,
        cancer_type="Melanoma",
        platform="RNA-seq",
        organ=Some("Skin"),
        has_response_labels=True,
    ),
    "VanAllen": DatasetSpec(
        id="VanAllen",
        title="Van Allen et al. 2015 Melanoma CTLA-4",
        modality=Modality.BULK_RNA,
        cancer_type="Melanoma",
        platform="RNA-seq",
        organ=Some("Skin"),
        has_response_labels=True,
    ),
    "EGAD00001006631": DatasetSpec(
        id="EGAD00001006631",
        title="Genentech IMvigor210 Align ReadsPerGene Counts",
        modality=Modality.BULK_RNA,
        cancer_type="Urothelial Bladder",
        platform="RNA-seq",
        organ=Some("Bladder"),
        has_response_labels=True,
    ),
}


def get_dataset_spec(dataset_id: str) -> Maybe[DatasetSpec]:
    """Retrieve the immutable dataset metadata specification by ID."""
    return Some(DATASET_REGISTRY[dataset_id]) if dataset_id in DATASET_REGISTRY else Nothing


def list_registered_datasets() -> tuple[DatasetSpec, ...]:
    """Return tuple of all registered dataset specifications."""
    return tuple(DATASET_REGISTRY.values())
