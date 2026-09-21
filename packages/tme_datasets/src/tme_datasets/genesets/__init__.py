"""Gene set collections, GMT parsers, overlap diagnostics, and signature scoring."""

from .collections import (
    AYERS_T_CELL_INFLAMED_GEP,
    TME_MAJOR_MARKERS,
    TME_SUBTYPE_MARKERS,
    get_bagaev_core_collection,
    get_tme_major_lineage_collection,
    get_tme_subtype_collection,
)
from .models import GeneSet, GeneSetCollection
from .overlap import compute_geneset_overlap
from .parser import export_gmt, parse_gmt
from .scoring import score_geneset_auc, score_geneset_zscore

__all__ = [
    "GeneSet",
    "GeneSetCollection",
    "AYERS_T_CELL_INFLAMED_GEP",
    "TME_MAJOR_MARKERS",
    "TME_SUBTYPE_MARKERS",
    "get_tme_major_lineage_collection",
    "get_tme_subtype_collection",
    "get_bagaev_core_collection",
    "parse_gmt",
    "export_gmt",
    "score_geneset_zscore",
    "score_geneset_auc",
    "compute_geneset_overlap",
]
