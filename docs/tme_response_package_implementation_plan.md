# Architecture & Implementation Plan: `tme_response` Package

## 1. Package Overview & Objectives

The `tme_response` package is a dedicated, functional-first Python library designed to predict and benchmark patient response to immune checkpoint therapies (ICT / ICB) by integrating bulk transcriptomic and genomic data.

It directly depends on the workspace package `tme_datasets` for access to standardized bulk RNA-seq cohorts, single-cell reference atlases (for pseudobulk simulation), curated gene sets, and clinical annotations.

```
                           UV Workspace Root: tme_analysis
                                         │
                 ┌───────────────────────┴───────────────────────┐
                 ▼                                               ▼
         packages/tme_datasets                         packages/tme_response
         ─────────────────────                         ─────────────────────
         • Raw cohort downloaders & caches             • Multi-omic ICB predictors (CYT, TIDE, IMPRES, GEP)
         • Bulk RNA expression matrices                • Genomic feature extractors (TMB, CNV, driver SNVs)
         • Clinical curation (RECIST, timepoints)      • DNA-RNA composite synergy models
         • Single-cell reference pools                 • Deconvolution & pseudobulk benchmark validation
         • Simulation (simulate_pseudobulk)            • Cross-cohort evaluation & Altair SVG reports
```

### Key Design Pillars
1. **Dependency on `tme_datasets`**: Uses `tme_datasets.load_dataset` as the authoritative single source of truth for expression matrices and clinical metadata.
2. **Multi-Omics Integration (RNA + DNA)**: Pairs transcriptomic inflamed/exhausted/excluded states with genomic determinants (Tumor Mutational Burden [TMB], frameshift indels, aneuploidy, and loss-of-heterozygosity [LOH] / driver mutations like *B2M*, *JAK1/2*, *STK11*, *KEAP1*).
3. **Strict Functional Programming (FP)**:
   * Pure mathematical functions mapping inputs to outputs without in-place mutation.
   * Monadic error handling using `Result` (`Success`/`Failure`) and `Maybe` (`Some`/`Nothing`) from `returns`.
   * Monadic pipeline composition using `@do`.
   * Immutable data models via Pydantic (`frozen=True`) or `NamedTuple`.
   * Declarative data processing via **Polars** (strictly avoiding pandas).
   * Publication-grade vector visualizations exported via **Altair** as **SVG**.
4. **Exhaustive Multi-Cohort Benchmarking**: Automated execution across all 9+ bulk ICB cohorts registered in `tme_datasets` with zero data leakage.

---

## 2. Directory Structure & Module Layout

```
packages/tme_response/
├── pyproject.toml
├── README.md
├── src/
│   └── tme_response/
│       ├── __init__.py
│       ├── types.py             # Enums: PredictorCategory, ModalityRequirement, ValidationMetric
│       ├── models.py            # Frozen Pydantic schemas: GenomicFeatures, ResponsePrediction, CohortMetrics
│       ├── data/
│       │   ├── __init__.py
│       │   ├── adapter.py       # Bridges tme_datasets AnnData -> MultiOmicCohort
│       │   ├── genomics.py      # MAF parser, TMB calculator, CNA/aneuploidy extraction
│       │   └── filters.py       # Baseline (Pre-treatment) filtering and RECIST binarization
│       ├── signatures/
│       │   ├── __init__.py
│       │   ├── cyt.py           # Rooney Cytolytic Activity (GZMA + PRF1)
│       │   ├── impres.py        # Auslander 15-pair non-parametric score
│       │   ├── gep.py           # Ayers 18-gene T-cell inflamed signature
│       │   ├── single_gene.py   # Litchfield CXCL9, CD8A, PDCD1, IFNG
│       │   ├── ipres.py         # Hugo 26-gene set innate resistance (EMT/angiogenesis)
│       │   └── scoring.py       # Unified signature scoring dispatcher
│       ├── models/
│       │   ├── __init__.py
│       │   ├── tide.py          # TIDE wrapper (zero-centering, dysfunction, exclusion)
│       │   ├── easier_bridge.py # Interop bridge for EaSIeR (PROGENy + DoRothEA)
│       │   └── iobr_bridge.py   # Interop bridge for IOBR 250+ signature suite
│       ├── deconvolution/
│       │   ├── __init__.py
│       │   ├── wrapper.py       # MCP-counter, EPIC, quanTIseq, CIBERSORTx
│       │   └── validation.py    # Pseudobulk benchmark against tme_datasets.simulation
│       ├── synergy/
│       │   ├── __init__.py
│       │   ├── dna_rna.py       # Composite scores: alpha*RNA + beta*log10(TMB)
│       │   └── gating.py        # Biological logic gating (e.g. CXCL9-high AND B2M-WT)
│       ├── evaluation/
│       │   ├── __init__.py
│       │   ├── metrics.py       # ROC-AUC, PR-AUC, Brier score, odds ratios
│       │   ├── survival.py      # Concordance Index, Cox Proportional Hazards for PFS/OS
│       │   └── benchmark.py     # Cross-cohort orchestrator across all tme_datasets
│       └── visualization/
│           ├── __init__.py
│           ├── roc_pr.py        # Altair ROC and PR curves (SVG export)
│           ├── benchmark_bar.py # Cross-cohort performance comparison bar charts
│           └── heatmap.py       # Signature collinearity and correlation matrix
└── tests/
    ├── test_signatures.py       # Unit tests for pure signature functions
    ├── test_genomics.py         # TMB and mutation parsing tests
    ├── test_synergy.py          # Composite model tests
    ├── test_adapter.py          # Ingestion tests with tme_datasets
    └── test_properties.py       # Hypothesis property-based tests
```

---

## 3. Detailed Architectural Components & Implementation Specifications

### Phase 1: Package Scaffolding & Dependencies (`pyproject.toml`)

The package will be declared within the uv workspace:

```toml
[project]
name = "tme_response"
version = "0.1.0"
description = "Multi-omic immune checkpoint response prediction and systematic benchmarking framework"
requires-python = ">=3.11"
dependencies = [
    "tme_datasets",              # Direct workspace dependency
    "anndata>=0.10.0",
    "polars>=0.20.0",
    "returns>=0.22.0",
    "pydantic>=2.5.0",
    "scipy>=1.11.0",
    "numpy>=1.24.0",
    "scikit-learn>=1.3.0",
    "altair>=5.2.0",
    "vl-convert-python>=1.2.0",
    "tidepy>=1.3.4",
]

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"
```

---

### Phase 2: Data & Genomic Ingestion Adapter (`tme_response.data`)

#### 2.1 MultiOmicCohort Data Structure
A frozen Pydantic model representing a synchronized cohort with aligned expression and genomic data:

```python
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe
import polars as pl

class MultiOmicCohort(BaseModel):
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    cohort_id: str
    cancer_type: str
    sample_ids: tuple[str, ...]
    expression_tpm: pl.DataFrame          # samples x genes
    response_binary: pl.DataFrame         # sample_id, response_binary (0/1), recist, timepoint
    tmb_scores: Maybe[pl.DataFrame]       # sample_id, tmb_per_mb, nonsynonymous_count
    driver_mutations: Maybe[pl.DataFrame] # sample_id, gene, variant_classification
    aneuploidy_scores: Maybe[pl.DataFrame]# sample_id, aneuploidy_score, arm_loss_count
```

#### 2.2 Genomic Extraction from `tme_datasets` Cache
While `tme_datasets` loads expression into `adata.X` and clinical data into `adata.obs`, the underlying downloaded archive (e.g. `data/raw/Liu-iAtlas/data_mutations.txt` or `data_cna.txt`) contains rich genomic tables:

```python
def extract_genomic_features(cohort_dir: Path, sample_ids: tuple[str, ...]) -> Result[GenomicTables, str]:
    """Pure parser extracting TMB and driver mutations from cBioPortal/MAF files."""
    # 1. Parse data_mutations.txt
    # Compute TMB = count of non-synonymous mutations / exome size (e.g. 38 Mb)
    # 2. Extract key resistance mutations (B2M, JAK1, JAK2, STK11, KEAP1, PTEN, EGFR)
    # 3. Parse data_cna.txt to compute 9p21/CDKN2A homozygous deletion and chromosome arm loss
```

#### 2.3 Baseline Pre-Treatment Filtering
Calls `tme_datasets.preprocessing.harmonize_obs_metadata` and filters strictly for:
* `biopsy_timepoint == "Pre"`
* `response_binary.is_not_null()`

---

### Phase 3: Pure Functional Signature Engine (`tme_response.signatures`)

All signature scoring functions are pure, non-mutating functions taking a `pl.DataFrame` or `MultiOmicCohort` and returning `Result[pl.DataFrame, str]`:

#### 3.1 Rooney CYT Activity (`signatures/cyt.py`)
$$\text{CYT} = \frac{\log_2(\text{GZMA} + 1) + \log_2(\text{PRF1} + 1)}{2}$$

#### 3.2 Auslander IMPRES (`signatures/impres.py`)
Evaluates 15 pairwise relations between immune checkpoint molecules.
Vectorized implementation using Polars expressions:
$$\text{IMPRES} = \sum_{k=1}^{15} \mathbb{I}(X_{A_k} > X_{B_k})$$

#### 3.3 Ayers Expanded T-Cell Inflamed GEP (`signatures/gep.py`)
Directly integrates with `tme_datasets.genesets.AYERS_T_CELL_INFLAMED_GEP` and computes standardized z-score or weighted mean.

#### 3.4 Litchfield Single-Gene Biomarkers (`signatures/single_gene.py`)
Extracts $\log_2(\text{TPM} + 1)$ for *CXCL9*, *CD8A*, *PDCD1*, and *IFNG*.

#### 3.5 Hugo IPRES (`signatures/ipres.py`)
Calculates single-sample enrichment across the 26 innate resistance gene sets (mesenchymal, angiogenesis, wound healing).

---

### Phase 4: TIDE & Systems Modeling Integrations (`tme_response.models`)

#### 4.1 TIDE Engine (`models/tide.py`)
* Computes cohort-wide zero-centering:
  $$\tilde{x}_{ij} = \log_2(x_{ij} + 1) - \mu_i$$
* Executes TIDE algorithm (dysfunction against CTL markers *CD8A, CD8B, GZMA, GZMB, PRF1*; exclusion against CAF, MDSC, M2 TAM profiles).
* Returns `TIDE_score`, `dysfunction_score`, `exclusion_score`.

#### 4.2 EaSIeR Bridge (`models/easier_bridge.py`)
* Exports expression to temporary Arrow/Parquet IPC buffer.
* Invokes `easier` via `rpy2` or headless CLI subprocess to extract:
  * PROGENy pathway scores (MAPK, NFkB, TGF-β, etc.).
  * DoRothEA transcription factor activities.
  * Predicted response probability.
* Reads back into Polars.

---

### Phase 5: Multi-Omic DNA + RNA Synergy Engine (`tme_response.synergy`)

A major finding of independent benchmarks (e.g. *Litchfield et al., Cell 2021*) is that RNA (inflamed TME) and DNA (tumor foreignness) are orthogonal, complementary determinants of response.

The `tme_response.synergy` module implements:

#### 5.1 Composite DNA-RNA Linear Index
$$S_{\text{composite}} = z(\text{RNA\_Score}) + \gamma \cdot z(\log_{10}(\text{TMB} + 1))$$
where $z(\cdot)$ denotes cohort-standardized z-scores, and $\gamma$ is a tunable or cross-validated weighting factor.

#### 5.2 Mechanistic Gating Logic (Biological Filters)
* **Antigen Presentation Gating**: Even in tumors with high *CXCL9* / GEP / CYT, inactivating mutations in *B2M* or *JAK1/2* abrogate response.
  $$\text{Score}_{\text{gated}} = \begin{cases} 
  \text{Score}_{\text{RNA}} & \text{if } B2M = \text{WT} \land JAK1 = \text{WT} \\
  \min(\text{Score}_{\text{RNA}}) & \text{if } B2M = \text{MUT} \lor JAK1 = \text{MUT}
  \end{cases}$$
* **Immune Exclusion by Copy Number Loss**: Combining high CAF / TIDE exclusion with 9p21 (*CDKN2A*) loss or high aneuploidy to identify deep immune-desert phenotypes.

---

### Phase 6: Deconvolution & Pseudobulk In Silico Validation (`tme_response.deconvolution`)

#### 6.1 Clinical Bulk Cohort Deconvolution
* Runs MCP-counter, EPIC, quanTIseq, and ESTIMATE across the bulk cohorts in `tme_datasets`.
* Extracts CD8+ T cell, Cytotoxic lymphocyte, Monocyte/Macrophage, and CAF fractions.

#### 6.2 Synthetic In Silico Validation Harness
* Ingests single-cell reference datasets from `tme_datasets` (`GSE120575` Sade-Feldman, `GSE115978` Jerby-Arnon, `Pelka_CRC`, etc.).
* Uses `tme_datasets.simulation.simulate_pseudobulk(adata, group_by="sample_id")` to generate simulated bulk samples where the true fractional cell-type counts are mathematically known.
* Measures deconvolution RMSE, Pearson $r$, and Lin's Concordance Correlation Coefficient ($CCC$) before applying to real clinical bulk cohorts.

---

### Phase 7: Systematic Evaluation & Reporting Engine (`tme_response.evaluation`)

#### 7.1 Cross-Cohort Validation Matrix
An automated benchmark orchestrator that runs all predictors across all registered bulk cohorts:

```python
COHORTS_TO_EVALUATE = (
    "Hugo-iAtlas",       # Melanoma (anti-PD-1)
    "Riaz-iAtlas",       # Melanoma (Nivolumab)
    "Liu-iAtlas",        # Melanoma (anti-PD-1)
    "Gide-iAtlas",       # Melanoma (anti-PD-1 +/- anti-CTLA-4)
    "Rosenberg-iAtlas",  # Bladder (Atezolizumab)
    "McDermott-iAtlas",  # RCC (Atezolizumab +/- Bevacizumab)
    "Padron-iAtlas",     # Pancreatic (anti-PD-1 + CD40)
    "Anders-iAtlas",     # TNBC (Atezolizumab)
    "VanAllen",          # Melanoma (anti-CTLA-4)
)
```

#### 7.2 Metric Computation
* **Discrimination**: ROC-AUC, Precision-Recall AUC (PR-AUC), Brier Score.
* **Survival Association**: Cox Proportional Hazards Hazard Ratio (HR) and log-rank p-value for PFS and OS.
* **Stability Index**: Coefficient of variation (CV) of AUC across independent cohorts to penalize overfitted signatures.

#### 7.3 Altair Vector Visualizations (SVG Export)
Strictly adheres to `.agents/rules/code-style-guide.md`:
* Generates ROC curves, PR curves, cross-cohort performance bar charts, and signature correlation heatmaps.
* Exports directly to **SVG** format using `vl-convert-python`.

---

## 4. Step-by-Step Execution Plan

| Step | Milestone | Deliverables / Files Created | Dependencies |
| :--- | :--- | :--- | :--- |
| **1** | **Scaffolding & Workspace Registration** | `packages/tme_response/pyproject.toml`, `src/tme_response/__init__.py`, `types.py`, `models.py`. Update root `pyproject.toml`. | None |
| **2** | **Data Adapter & Genomic Extraction** | `data/adapter.py`, `data/genomics.py`, `data/filters.py`. Ingests `tme_datasets.load_dataset`. | Step 1 |
| **3** | **Pure Transcriptomic Signatures** | `signatures/cyt.py`, `impres.py`, `gep.py`, `single_gene.py`, `ipres.py`. Unit tests in `tests/test_signatures.py`. | Step 2 |
| **4** | **TIDE & Systems Models** | `models/tide.py` (with cohort zero-centering), `models/easier_bridge.py`. | Step 2 |
| **5** | **DNA-RNA Synergy Models** | `synergy/dna_rna.py` (composite scoring), `synergy/gating.py` (*B2M* / *JAK1* gating). | Steps 2 & 3 |
| **6** | **Deconvolution & In Silico Validation** | `deconvolution/wrapper.py`, `deconvolution/validation.py` using `tme_datasets.simulation.simulate_pseudobulk`. | Steps 1 & 2 |
| **7** | **Evaluation & Benchmark Suite** | `evaluation/metrics.py`, `evaluation/benchmark.py`, `visualization/roc_pr.py`. Automated runner across 9 cohorts. | Steps 3–6 |
| **8** | **Property-Based Testing & Verification** | Hypothesis tests (`tests/test_properties.py`), `mypy --strict` compliance. | Steps 1–7 |

---

## 5. Verification & Testing Strategy

1. **Unit & Functional Tests (`pytest`)**:
   * Test signature calculations against analytical toy matrices.
   * Verify that IMPRES produces integer scores in $[0, 15]$.
   * Verify that CYT equals the geometric mean of *GZMA* and *PRF1*.
2. **Property-Based Testing (`hypothesis`)**:
   * Invariant: Scaling all gene expressions by a positive constant must not alter rank-based scores (IMPRES rank invariance).
   * Invariant: Adding zero-expression genes must not crash or distort signature scoring.
3. **Static Analysis & Purity**:
   * Strict typing via `mypy --strict`.
   * Enforce no `None` / `Optional` (using `Maybe[T]`).
   * Enforce no exception-based control flow (using `Result[T, str]`).
4. **End-to-End Benchmark Execution**:
   * Run benchmark across `"Hugo-iAtlas"` and `"Rosenberg-iAtlas"`.
   * Verify that all ROC-AUC and PR-AUC scores are produced without errors and output to SVG via Altair.
