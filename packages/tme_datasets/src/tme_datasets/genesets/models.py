"""Data models for gene sets, collections, and signature scoring."""

from __future__ import annotations

from typing import Mapping
from pydantic import BaseModel, ConfigDict
import polars as pl


class GeneSet(BaseModel):
    """A biological gene set or signature."""

    model_config = ConfigDict(frozen=True)

    id: str
    name: str
    description: str = ""
    genes: tuple[str, ...]
    organism: str = "human"


class GeneSetCollection(BaseModel):
    """A curated collection of gene sets (e.g. Bagaev MFP, MSigDB, ImmunoCompass)."""

    model_config = ConfigDict(frozen=True)

    id: str
    name: str
    description: str = ""
    version: str = "1.0"
    gene_sets: Mapping[str, GeneSet]

    def get_gene_set(self, key: str) -> GeneSet | None:
        return self.gene_sets.get(key)

    @property
    def all_genes(self) -> tuple[str, ...]:
        unique = sorted({g for gs in self.gene_sets.values() for g in gs.genes})
        return tuple(unique)
