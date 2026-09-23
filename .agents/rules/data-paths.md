# Mandatory Data Paths Resolution Standard

This document outlines the strict rules for resolving file paths, storing datasets, and caching in the `tme_analysis` repository. All AI agents working in this codebase must strictly adhere to these standards.

## 1. Zero Path Hardcoding Policy
- **Strictly Banned**: AI agents MUST NEVER hardcode path strings (such as `"data/preprocessed/..."`, `"data/raw/..."`, or `"scratch/..."`) in Python modules, Jupyter notebooks, scripts, or workflows.
- **Mandatory Query API**: All paths must be dynamically resolved through the centralized `tme_datasets.paths` module.

```python
# FORBIDDEN (DO NOT WRITE THIS):
h5ad_path = Path("data/preprocessed/GSE120575.h5ad")
raw_dir = Path("data/raw/GSE120575")

# REQUIRED (ALWAYS USE THIS):
from tme_datasets.paths import (
    get_data_paths,
    get_preprocessed_h5ad_path,
    get_raw_dataset_dir,
    get_scratch_dataset_dir,
    find_dataset_h5ad,
)

h5ad_path = get_preprocessed_h5ad_path("GSE120575")
raw_dir = get_raw_dataset_dir("GSE120575")
```

## 2. Deterministic Repository Layout
All paths are backed by `config/data_paths.toml` and adhere to this hierarchy:
```text
/Users/halu/Code/tme_analysis/
├── config/
│   └── data_paths.toml                    # Master data configuration (paths, templates, overrides)
├── data/
│   ├── ensembl/                           # PYENSEMBL LOCAL DATA & GTF DATABASES
│   │   ├── homo_sapiens/                  # Downloaded GTFs and SQLite indexes
│   │   └── gene_mapping_cache_release_111.parquet # Fast symbol-to-ENSG persistent cache
│   ├── raw/
│   │   ├── <dataset_id>/                  # Isolated directory for raw downloaded archives (tar.gz, txt.gz, csv.gz)
│   ├── preprocessed/
│   │   ├── <dataset_id>.h5ad              # CANONICAL PREPROCESSED ANNDATA OBJECTS (<0.5s loading)
│   └── reference.h5ad                     # Combined multi-cohort reference atlas
├── dataset_papers/                        # Read-only published paper H5AD files
└── scratch/                               # Ephemeral intermediate shards, tests, logs
```

## 3. Fast-Loading & H5AD Caching Invariant
- **Always Prioritize H5AD**: `load_dataset("<dataset_id>")` checks `find_dataset_h5ad("<dataset_id>")` first. When an `.h5ad` file exists, it MUST load directly in `<0.5s` and NEVER access or decompress raw `.txt` / `.tsv` files.
- **Automatic Caching**: When ingesting new or raw cohorts, once parsed into an `AnnData` object, it MUST immediately be written to `get_preprocessed_h5ad_path("<dataset_id>")` so that all subsequent invocations load from H5AD.
- **Cache Invalidation**:
  - `load_dataset(id, force_recompute=True)`: Re-parses raw text data and overwrites the `.h5ad` cache without re-downloading existing raw files.
  - `load_dataset(id, force_download=True)`: Re-downloads raw vendor archives from network before parsing.

## 4. Jupyter Notebook Compatibility
- When writing scripts or notebook cells, use `load_dataset` which flushes logs directly to `sys.stdout` so progress is visible in real-time without buffering delay or stderr warning boxes.
