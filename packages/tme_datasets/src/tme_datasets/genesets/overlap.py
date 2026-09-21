"""Pairwise overlap and Jaccard similarity between gene sets."""

from __future__ import annotations

import polars as pl
from returns.result import Failure, Result, Success

from .models import GeneSetCollection


def compute_geneset_overlap(
    collection_a: GeneSetCollection,
    collection_b: GeneSetCollection,
) -> Result[pl.DataFrame, str]:
    """Calculate pairwise Jaccard index and common gene counts between two collections."""
    try:
        rows = []
        for id_a, gs_a in collection_a.gene_sets.items():
            set_a = set(gs_a.genes)
            for id_b, gs_b in collection_b.gene_sets.items():
                set_b = set(gs_b.genes)
                intersection = len(set_a & set_b)
                union = len(set_a | set_b)
                jaccard = intersection / union if union > 0 else 0.0

                rows.append({
                    "signature_a": id_a,
                    "signature_b": id_b,
                    "n_genes_a": len(set_a),
                    "n_genes_b": len(set_b),
                    "n_shared_genes": intersection,
                    "jaccard_similarity": round(jaccard, 4),
                })

        return Success(pl.DataFrame(rows))
    except Exception as exc:
        return Failure(f"Failed to compute gene set overlap: {exc}")
