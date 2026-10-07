"""Standardized clinical single-cell cohort definitions and metadata mappings.
"""

from __future__ import annotations

# Standardized Clinical Single-Cell Cohort Layout (3x3 grid, panels a-i)
COHORT_PANELS: tuple[dict[str, str], ...] = (
    {"tag": "a", "acc": "GSE120575", "indication": "Melanoma", "tech": "Smart-seq2", "patients": "14 R / 25 NR"},
    {"tag": "b", "acc": "CELLxGENE_7b20c613", "indication": "Melanoma", "tech": "10x Chromium 3'", "patients": "41 R / 82 NR"},
    {"tag": "c", "acc": "CELLxGENE_05a8c945", "indication": "Colorectal", "tech": "10x Chromium 3'/5'", "patients": "48 R / 29 NR"},
    {"tag": "d", "acc": "CELLxGENE_6f9de485", "indication": "Breast", "tech": "10x Chromium 3'", "patients": "45 R / 37 NR"},
    {"tag": "e", "acc": "GSE207422", "indication": "NSCLC", "tech": "High-Throughput", "patients": "4 R / 10 NR"},
    {"tag": "f", "acc": "GSE243013", "indication": "NSCLC", "tech": "High-Throughput", "patients": "130 R / 112 NR"},
    {"tag": "g", "acc": "GSE233203", "indication": "NSCLC", "tech": "10x Chromium 5'", "patients": "3 R / 4 NR"},
    {"tag": "h", "acc": "GSE200996", "indication": "HNSCC", "tech": "High-Throughput", "patients": "4 R / 8 NR"},
    {"tag": "i", "acc": "GSE316195", "indication": "PDAC", "tech": "Single-Nucleus", "patients": "15 R / 2 NR"},
)

COHORT_METADATA: dict[str, dict[str, str]] = {
    p["acc"]: {"indication": p["indication"], "tech": p["tech"], "patients": p["patients"]}
    for p in COHORT_PANELS
}
