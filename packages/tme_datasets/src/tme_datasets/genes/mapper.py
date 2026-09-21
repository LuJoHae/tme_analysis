"""Gene identifier normalization and mapping engine."""

from __future__ import annotations

import re
from typing import Mapping
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success

from ..types import GeneIDType

# Common TME checkpoint and lineage gene mapping table for instant offline resolution
OFFLINE_ENSEMBL_TO_HUGO: Mapping[str, str] = {
    "ENSG00000012048": "BRCA1",
    "ENSG00000139618": "BRCA2",
    "ENSG00000188389": "PDCD1",      # PD-1
    "ENSG00000120217": "CD274",      # PD-L1
    "ENSG00000163599": "CTLA4",      # CTLA-4
    "ENSG00000153563": "CD8A",
    "ENSG00000172116": "CD8B",
    "ENSG00000010610": "CD4",
    "ENSG00000198851": "CD3E",
    "ENSG00000167286": "CD3D",
    "ENSG00000011465": "CD68",
    "ENSG00000170458": "CD14",
    "ENSG00000105383": "CD19",
    "ENSG00000156738": "MS4A1",      # CD20
    "ENSG00000049768": "FOXP3",
    "ENSG00000026025": "VIM",
    "ENSG00000119888": "EPCAM",
    "ENSG00000115085": "ZAP70",
    "ENSG00000135046": "HAVCR2",     # TIM-3
    "ENSG00000089692": "LAG3",
    "ENSG00000181847": "TIGIT",
    "ENSG00000100385": "IL2RA",
    "ENSG00000107742": "COL1A1",
    "ENSG00000164692": "COL1A2",
    "ENSG00000105329": "PECAM1",     # CD31
    "ENSG00000110799": "VWF",
    "ENSG00000100097": "LGALS9",
    "ENSG00000163600": "ICOS",
    "ENSG00000115902": "NKG7",
    "ENSG00000149294": "NCAM1",      # CD56
    "ENSG00000180644": "PRF1",
    "ENSG00000100453": "GZMB",
    "ENSG00000113088": "GZMK",
    "ENSG00000105374": "NKG2A",
    "ENSG00000111537": "IFNG",
}

OFFLINE_HUGO_TO_ENSEMBL: Mapping[str, str] = {
    v: k for k, v in OFFLINE_ENSEMBL_TO_HUGO.items()
}


def strip_gene_version(identifier: str) -> str:
    """Strip transcript or gene version decimal suffix (e.g. ENSG00000139618.15 -> ENSG00000139618)."""
    return identifier.split(".")[0] if identifier.startswith("ENS") and "." in identifier else identifier


def detect_gene_id_type(sample_identifiers: tuple[str, ...]) -> GeneIDType:
    """Infer the gene nomenclature format from a representative sample of identifiers."""
    if not sample_identifiers:
        return GeneIDType.HUGO_SYMBOL

    clean_samples = tuple(strip_gene_version(s) for s in sample_identifiers[:50])
    ensembl_matches = sum(1 for s in clean_samples if re.match(r"^ENSG\d{11}$", s))
    entrez_matches = sum(1 for s in clean_samples if s.isdigit())

    if ensembl_matches / len(clean_samples) > 0.4:
        return GeneIDType.ENSEMBL_ID
    if entrez_matches / len(clean_samples) > 0.4:
        return GeneIDType.ENTREZ_ID
    return GeneIDType.HUGO_SYMBOL


def map_gene_identifier(
    identifier: str,
    target_type: GeneIDType,
    strip_suffix: bool = True,
) -> Maybe[str]:
    """Map a single gene identifier to the target type using cached mappings or heuristics."""
    clean_id = strip_gene_version(identifier) if strip_suffix else identifier

    match target_type:
        case GeneIDType.HUGO_SYMBOL:
            if not clean_id.startswith("ENSG"):
                return Some(clean_id)
            return (
                Some(OFFLINE_ENSEMBL_TO_HUGO[clean_id])
                if clean_id in OFFLINE_ENSEMBL_TO_HUGO
                else Nothing
            )
        case GeneIDType.ENSEMBL_ID:
            if clean_id.startswith("ENSG"):
                return Some(clean_id)
            return (
                Some(OFFLINE_HUGO_TO_ENSEMBL[clean_id.upper()])
                if clean_id.upper() in OFFLINE_HUGO_TO_ENSEMBL
                else Nothing
            )
        case GeneIDType.ENTREZ_ID | GeneIDType.AUTO_DETECT:
            return Some(clean_id)
