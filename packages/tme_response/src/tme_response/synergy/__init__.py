"""Multi-omic synergy models integrating genomics and transcriptomics."""

from .dna_rna import compute_dna_rna_composite
from .gating import apply_antigen_presentation_gating

__all__ = [
    "compute_dna_rna_composite",
    "apply_antigen_presentation_gating",
]
