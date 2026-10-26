"""Standardized clinical single-cell cohort definitions and metadata mappings.
"""

from __future__ import annotations

# Standardized Clinical Single-Cell Cohort Layout (3x3 grid, panels a-i)
COHORT_PANELS: tuple[dict[str, str], ...] = (
    {"tag": "a", "acc": "GSE120575", "indication": "Melanoma", "author": "Sade-Feldman et al.", "year": "2018", "journal": "Cell", "tech": "Smart-seq2", "patients": "14 R / 25 NR"},
    {"tag": "b", "acc": "CELLxGENE_7b20c613", "indication": "Pan-Cancer", "author": "Gondal et al.", "year": "2025", "journal": "Sci Data", "tech": "10x Chromium 3'", "patients": "41 R / 82 NR"},
    {"tag": "c", "acc": "CELLxGENE_05a8c945", "indication": "Colorectal", "author": "Marteau et al.", "year": "2026", "journal": "Cancer Cell", "tech": "10x Chromium 3'/5'", "patients": "48 R / 29 NR"},
    {"tag": "d", "acc": "CELLxGENE_6f9de485", "indication": "Breast", "author": "Yan et al.", "year": "2026", "journal": "Nature", "tech": "10x Chromium 3'", "patients": "45 R / 37 NR"},
    {"tag": "e", "acc": "GSE207422", "indication": "NSCLC", "author": "Hu et al.", "year": "2023", "journal": "Genome Med", "tech": "High-Throughput", "patients": "4 R / 10 NR"},
    {"tag": "f", "acc": "GSE243013", "indication": "NSCLC", "author": "Liu et al.", "year": "2025", "journal": "Cell", "tech": "High-Throughput", "patients": "130 R / 112 NR"},
    {"tag": "g", "acc": "GSE233203", "indication": "NSCLC", "author": "Hong et al.", "year": "2025", "journal": "MedComm", "tech": "10x Chromium 5'", "patients": "3 R / 4 NR"},
    {"tag": "h", "acc": "GSE200996", "indication": "HNSCC", "author": "Luoma et al.", "year": "2022", "journal": "Cell", "tech": "High-Throughput", "patients": "4 R / 8 NR"},
    {"tag": "i", "acc": "GSE316195", "indication": "PDAC", "author": "Raufi et al.", "year": "2026", "journal": "Nat Commun", "tech": "Single-Nucleus", "patients": "15 R / 2 NR"},
)

COHORT_METADATA: dict[str, dict[str, str]] = {
    p["acc"]: {
        "indication": p["indication"],
        "author": p["author"],
        "year": p["year"],
        "journal": p["journal"],
        "tech": p["tech"],
        "patients": p["patients"],
    }
    for p in COHORT_PANELS
}


def get_single_cell_cohort_display_label(
    accession: str,
    format_style: str = "cancer_author_year",
) -> str:
    """Returns publication-grade display label for single-cell response cohorts.

    Format options:
    - 'cancer_author_year': 'Melanoma · Sade-Feldman et al. (2018)' (default)
    - 'author_year_cancer': 'Sade-Feldman et al. (2018) [Melanoma]'
    - 'comprehensive': 'Melanoma: Sade-Feldman et al. (2018) [GSE120575]'
    """
    meta = COHORT_METADATA.get(accession)
    if meta is None:
        return accession

    cancer = meta["indication"]
    author = meta["author"]
    year = meta["year"]

    match format_style:
        case "cancer_author_year_3line":
            return f"{cancer} · {author} · ({year})"
        case "cancer_author_year":
            return f"{cancer} · {author} ({year})"
        case "author_year_cancer":
            return f"{author} ({year}) [{cancer}]"
        case "comprehensive":
            return f"{cancer}: {author} ({year}) [{accession}]"
        case _:
            return f"{cancer} · {author} · ({year})"

