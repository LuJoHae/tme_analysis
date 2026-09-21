"""Confounding gene family filtering for reference and deconvolution signatures."""

from __future__ import annotations

import re
import anndata as ad
from returns.result import Failure, Result, Success

CONFOUNDING_PATTERNS: tuple[str, ...] = (
    r"^RP[SL]\d+",           # Ribosomal proteins
    r"^MT-",                 # Mitochondrial
    r"^IG[HKL]",             # Immunoglobulin chains
    r"^TR[ABGD][CVJD]",      # T-cell receptor variable/constant chains
    r"^RP11-|^AC\d{6}-",     # Uncharacterized clones/pseudogenes
    r"^RNA5S|^RNU",          # Small RNAs
    r"^LINC\d+",             # Long intergenic non-coding
    r"^(MLANA|PMEL|TYR|DCT|MITF|S100B|MAGEA\d+)",  # Melanoma / tumor lineage markers
)


def filter_confounding_genes(
    adata: ad.AnnData,
    patterns: tuple[str, ...] = CONFOUNDING_PATTERNS,
) -> Result[ad.AnnData, str]:
    """Filter out confounding, technical, and lineage-specific artifact genes from AnnData."""
    try:
        combined_regex = re.compile("|".join(patterns), flags=re.IGNORECASE)
        var_names = [str(g) for g in adata.var_names]
        keep_mask = [not bool(combined_regex.search(g)) for g in var_names]

        filtered = adata[:, keep_mask].copy()
        return Success(filtered)
    except Exception as exc:
        return Failure(f"Failed to filter confounding genes: {exc}")
