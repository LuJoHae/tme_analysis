"""Declarative GMT and tabular parser for gene set collections."""

from __future__ import annotations

from pathlib import Path
from returns.result import Failure, Result, Success

from .models import GeneSet, GeneSetCollection


def parse_gmt(file_path: Path, collection_id: str = "custom") -> Result[GeneSetCollection, str]:
    """Parse a GMT (Gene Matrix Transposed) file into an immutable GeneSetCollection."""
    if not file_path.is_file():
        return Failure(f"GMT file not found: {file_path}")

    try:
        gene_sets = {}
        with open(file_path, "r", encoding="utf-8") as fh:
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    set_id = parts[0].strip()
                    desc = parts[1].strip()
                    genes = tuple(filter(None, (g.strip() for g in parts[2:])))
                    if set_id and genes:
                        gene_sets[set_id] = GeneSet(
                            id=set_id,
                            name=set_id,
                            description=desc,
                            genes=genes,
                        )

        return Success(
            GeneSetCollection(
                id=collection_id,
                name=file_path.stem,
                gene_sets=gene_sets,
            )
        )
    except Exception as exc:
        return Failure(f"Failed to parse GMT file {file_path}: {exc}")


def export_gmt(collection: GeneSetCollection, out_path: Path) -> Result[Path, str]:
    """Export a GeneSetCollection to standard GMT format."""
    try:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w", encoding="utf-8") as fh:
            for gs in collection.gene_sets.values():
                genes_str = "\t".join(gs.genes)
                fh.write(f"{gs.id}\t{gs.description}\t{genes_str}\n")
        return Success(out_path)
    except Exception as exc:
        return Failure(f"Failed to export GMT to {out_path}: {exc}")
