# Workspace Package Map

This repository is organized as a monorepo managed by `uv`. All core domain logic, algorithms, models, and data accessors live inside modular Python packages under `packages/*` and are installed in the workspace environment as editable packages.

When writing or modifying code, agents must reference and reuse these packages rather than implementing local routines in `scripts/` or `workflow/`.

---

## Monorepo Architecture & Dependency Hierarchy

To maintain clean separation of concerns, prevent circular imports, and ensure modularity, packages follow a strict layered hierarchy:

1. **Layer 0: Foundational Data Layer (Strict Zero Internal Dependencies)**:
   - `tme_datasets`: The central data framework, storage, and registry. **Must NEVER depend on or import any other package in this repository.** All its dependencies must be standard external PyPI packages declared in its `pyproject.toml`.
2. **Layer 1: Specialized Biological & Algorithmic Packages**:
   - `gene_utils`, `bayesprism`, `instaprism`, `milopy`, `references`, `selective_inference`, `single_cell_datasets`, `single_cell_immuno_datasets`, `ici_datasets`, `genentech_datasets`.
   - May import from Layer 0 (`tme_datasets`), but must never introduce circular cross-dependencies with each other.
3. **Layer 2: Pipelines, High-Level Models & Workflows**:
   - `ml_pipelines`, `scripts/`, `workflow/`.
   - May import and compose from Layer 0 and Layer 1.

---

## Workspace Packages Index

### 1. `gene_utils`
- **Location**: `packages/gene_utils`
- **Import**: `import gene_utils`
- **Domain**: Gene identifiers, mapping, genomic coordinates, and annotations.
- **Key Responsibilities**:
  - Ensembl ID to gene symbol conversions (and vice versa).
  - HGNC symbol normalization and alias resolution.
  - Coordinate extraction, GTF/GFF parsing utilities, and genomic lookups.

### 2. `tme_datasets`
- **Location**: `packages/tme_datasets`
- **Import**: `import tme_datasets`
- **Domain**: Central tumor microenvironment (TME) datasets, loaders, and preprocessing.
- **Key Responsibilities**:
  - Unified registry for bulk, single-cell, and spatial transcriptomics datasets.
  - Standardized dataset loaders, caching, and serialization (Parquet, AnnData, HDF5).
  - Cohort harmonization, sample metadata management, and cross-cohort merging.
  - Perturbation frameworks, gene sets, and signature scoring.

### 3. `single_cell_datasets`
- **Location**: `packages/single_cell_datasets`
- **Import**: `import single_cell_datasets`
- **Domain**: Single-cell transcriptomics dataset loaders and reference atlases.
- **Key Responsibilities**:
  - Reference single-cell RNA-seq datasets and formatting.
  - AnnData container creation and standardization for single-cell cohorts.

### 4. `single_cell_immuno_datasets`
- **Location**: `packages/single_cell_immuno_datasets`
- **Import**: `import single_cell_immuno_datasets`
- **Domain**: Immuno-oncology single-cell cohorts and differential abundance.
- **Key Responsibilities**:
  - Cohort loaders (e.g. GSE120575, Gondal 2025).
  - Automated downloading and caching pipelines for single-cell ICB cohorts.
  - Preprocessing pipelines specifically tailored for Milo differential abundance analysis.

### 5. `ici_datasets`
- **Location**: `packages/ici_datasets`
- **Import**: `import ici_datasets`
- **Domain**: Immune Checkpoint Inhibitor / Blockade (ICI/ICB) bulk clinical cohorts.
- **Key Responsibilities**:
  - Bagaev et al. dataset loaders and subtype definitions.
  - cBioPortal clinical cohort fetching, curation, and standardization.
  - Response metadata (RECIST, progression-free survival, overall survival).

### 6. `genentech_datasets`
- **Location**: `packages/genentech_datasets`
- **Import**: `import genentech_datasets`
- **Domain**: Genentech-specific trial data and molecular profiles.
- **Key Responsibilities**:
  - Specialized clinical trial cohorts (e.g. IMvigor210, etc.) and accompanying expression/clinical data.

### 7. `ml_pipelines`
- **Location**: `packages/ml_pipelines`
- **Import**: `import ml_pipelines`
- **Domain**: Machine learning, feature selection, evaluation, and mutation profiling.
- **Key Responsibilities**:
  - Random Forest, logistic regression, and cross-validation pipelines.
  - Feature selection, ranking, and hyperparameter tuning (Optuna integration).
  - Tumor Mutation Burden (TMB), Variant Allele Fraction (VAF), and mutation matrix processing.
  - TCGA background distributions and comparative normalization.

### 8. `references`
- **Location**: `packages/references`
- **Import**: `import references`
- **Domain**: Canonical biological references, gene signatures, and marker catalogs.
- **Key Responsibilities**:
  - Cell-type marker gene lists (LM22, TME subtypes, Bagaev signatures).
  - Standardized reference signatures for deconvolution and pathway analysis.

### 9. `bayesprism`
- **Location**: `packages/bayesprism`
- **Import**: `import bayesprism`
- **Domain**: Bayesian cell type and gene expression deconvolution.
- **Key Responsibilities**:
  - Functional PyTorch implementation of BayesPrism.
  - Reference profile derivation, Gibbs sampling, optimization, and cell fraction posterior estimation.

### 10. `instaprism`
- **Location**: `packages/instaprism`
- **Import**: `import instaprism`
- **Domain**: High-speed linear and penalization deconvolution algorithms.
- **Key Responsibilities**:
  - Rapid cell-type fraction estimation from bulk expression profiles.

### 11. `milopy`
- **Location**: `packages/milopy`
- **Import**: `import milopy`
- **Domain**: Differential cell abundance testing on single-cell k-NN graphs.
- **Key Responsibilities**:
  - Native Python implementation of miloR.
  - k-NN graph neighborhood building, cell counting, GLM fitting, and spatial FDR correction.

### 12. `selective_inference`
- **Location**: `packages/selective_inference`
- **Import**: `import selective_inference`
- **Domain**: Selective statistical inference and post-selection hypothesis testing.
- **Key Responsibilities**:
  - Truncated Gaussian/chi-square conditioning, clustering-adjusted p-values, and simulation routines.
