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
    get_manual_download_dir,
    get_scratch_dataset_dir,
    find_dataset_h5ad,
)

h5ad_path = get_preprocessed_h5ad_path("GSE120575")
raw_dir = get_raw_dataset_dir("GSE120575")
manual_dir = get_manual_download_dir()
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
│   ├── manual_download/                   # MANUALLY DOWNLOADED & CONTROLLED DATASETS (EGAD, Maynard, Paper H5ADs)
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

## 5. Manual Downloads, Third-Party Isolation, & Ingestion Standards

### A. Zero Fallback on External Wrapper Packages
- **Never Fall Back on Wrapper Packages**: AI agents MUST NOT fall back on external data wrapper packages or foreign repositories to fetch core study cohorts.
- **Controlled & Manual Cohorts**: Datasets that cannot be fetched via standard automated HTTP/FTP/API endpoints (e.g., controlled-access EGA data, proprietary cohorts, or datasets requiring portal login) MUST reside in `data/manual_download/` (resolved via `get_manual_download_dir()`).
- **Informative Failure Monad**: When a dataset file in `data/manual_download/` is missing, loaders must return an informative `Failure` specifying:
  1. The exact missing file path;
  2. Clear instructions or commands to acquire or stage it (e.g. `python scripts/setup_manual_downloads.py`).

### B. Reusable Package Ingestion Functions (No Standalone Scripts)
- All downloaders (e.g., Google Drive, GEO fetchers), format converters, and AnnData builders MUST be implemented as reusable, tested functions inside `packages/tme_datasets/src/tme_datasets/providers/` and `tme_datasets/src/tme_datasets/download/`.
- They must be exported at package top-level and registered in `tme_datasets.load_dataset(id, auto_download=True)` and `registry.DATASET_REGISTRY`. Never write one-off conversion scripts in the root or `scratch/`.

### C. AnnData H5AD Serialization Guard (Python 3.14 / h5py)
- When writing `AnnData` to `.h5ad` (`adata.write_h5ad()`), metadata DataFrames (`adata.obs`, `adata.var`) must be sanitized:
  1. Drop or rename empty string column headers (`[c for c in df.columns if not c or c == ""]`), as `h5py` crashes with `TypeError` when serializing empty key strings on Python 3.14.
  2. Sanitize object columns to ensure string serialization without dangling `None` values (`df[col].fillna("").astype(str)`).

### D. Streaming Sparse Ingestion for Large Expression Matrices
- When parsing multi-gigabyte text matrices (`.csv`, `.tsv`), DO NOT use `pandas.read_csv()` to load a dense matrix into memory.
- Stream line-by-line using vectorized numpy routines (`np.fromstring(line, sep=",", dtype=np.float32)`, extracting non-zeros via `np.flatnonzero()`), and construct a `scipy.sparse.csr_matrix` directly to prevent memory exhaustion.
