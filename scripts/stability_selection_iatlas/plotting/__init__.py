"""Publication-grade plotting modules for iAtlas Stability Selection analysis.

Adheres strictly to the Nature Methods minimalist wireframe standard,
Okabe-Ito colorblind-safe palettes, and native Inkscape layer architecture.
"""

from __future__ import annotations

from .combine_all_stability_paths import (
    PANEL_CONFIGS,
    PanelConfig,
    assemble_composite_svg,
)
from .plot_cohort_gene_overlap import (
    compute_upset_intersections,
    extract_cohort_models,
    extract_matrix_cells,
    render_full_svg as render_cohort_overlap_svg,
)
from .plot_fitter_gene_overlap import (
    extract_fitter_hits,
    load_benchmark_data,
    render_full_svg as render_fitter_overlap_svg,
)
from .plot_stability_selection_iatlas import (
    build_cohort_paths_chart,
    export_chart_to_svg,
)

__all__ = [
    "PANEL_CONFIGS",
    "PanelConfig",
    "assemble_composite_svg",
    "build_cohort_paths_chart",
    "compute_upset_intersections",
    "export_chart_to_svg",
    "extract_cohort_models",
    "extract_fitter_hits",
    "extract_matrix_cells",
    "load_benchmark_data",
    "render_cohort_overlap_svg",
    "render_fitter_overlap_svg",
]
