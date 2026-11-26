"""In-silico pseudobulk mixture simulation and empirical sample/patient pseudobulk aggregation."""

from .empirical_pseudobulk import (
    aggregate_empirical_pseudobulk,
    cache_empirical_pseudobulk,
    compute_empirical_cell_fractions,
    empirical_pseudobulk_to_anndata,
    load_empirical_pseudobulk,
    resolve_adata_gene_symbols,
    resolve_sample_column_name,
)
from .sc_pseudobulk_stream import (
    DeduplicateMode,
    ScaleMode,
    TumorPseudobulkResult,
    TumorPseudobulkSpec,
    aggregate_gzip_h5ad_tumor_pseudobulk,
    aggregate_multi_file_tumor_pseudobulk,
    build_tumor_cell_mask,
    classify_broad_lineage,
    deduplicate_patient_profiles,
    extract_h5ad_var_symbols,
    load_ensembl_symbol_map,
    stream_cohort_tumor_pseudobulk,
    stream_tumor_pseudobulk_h5ad,
)
from .simulator import simulate_pseudobulk

__all__ = [
    "DeduplicateMode",
    "ScaleMode",
    "TumorPseudobulkResult",
    "TumorPseudobulkSpec",
    "aggregate_empirical_pseudobulk",
    "aggregate_gzip_h5ad_tumor_pseudobulk",
    "aggregate_multi_file_tumor_pseudobulk",
    "build_tumor_cell_mask",
    "cache_empirical_pseudobulk",
    "classify_broad_lineage",
    "compute_empirical_cell_fractions",
    "deduplicate_patient_profiles",
    "empirical_pseudobulk_to_anndata",
    "extract_h5ad_var_symbols",
    "load_empirical_pseudobulk",
    "load_ensembl_symbol_map",
    "resolve_adata_gene_symbols",
    "resolve_sample_column_name",
    "simulate_pseudobulk",
    "stream_cohort_tumor_pseudobulk",
    "stream_tumor_pseudobulk_h5ad",
]
