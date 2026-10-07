# `tme_datasets` Comprehensive Interactive Tutorial

Welcome to the comprehensive tutorial for **`tme_datasets`**, the unified data access, downloading, preprocessing, harmonization, and perturbation framework for tumor microenvironment (TME) analysis.

This tutorial guides you through all key features step-by-step with practical, runnable examples.

---

## Quick Run

You can execute the complete end-to-end tutorial script right now from your terminal:

```bash
uv run python packages/tme_datasets/examples/tutorial_explore_all_features.py
```

---

## Tutorial Table of Contents

1. [Dataset Registry & Metadata Inspection](#1-dataset-registry--metadata-inspection)
2. [Loading Individual Datasets](#2-loading-individual-datasets)
3. [Gene Set Collections & Signature Scoring](#3-gene-set-collections--signature-scoring)
4. [Subsampling, Supersampling & Perturbations](#4-subsampling-supersampling--perturbations)
5. [In-Silico Pseudobulk Deconvolution Simulator](#5-in-silico-pseudobulk-deconvolution-simulator)
6. [Multi-Dataset Harmonization & Integration Metrics](#6-multi-dataset-harmonization--integration-metrics)
7. [PyTorch Dataset & DataLoader Bridge](#7-pytorch-dataset--dataloader-bridge)
8. [Gene Identifier Unification](#8-gene-identifier-unification)
9. [Out-of-Core Storage & Cryptographic Verification](#9-out-of-core-storage--cryptographic-verification)
10. [End-to-End Single-Cell Downstream Workflow (Case Study: GSE120575)](#10-end-to-end-single-cell-downstream-workflow-case-study-gse120575)

---

## 1. Dataset Registry & Metadata Inspection

All single-cell references, bulk iAtlas cohorts, and direct publication datasets are declaratively indexed in an immutable registry. Calling `list_registered_datasets()` returns a **`RegisteredDatasets`** collection that inherits from `tuple[DatasetSpec, ...]`, while providing powerful list accessors, dual indexing, and seamless export to **Polars** DataFrames.

```python
import polars as pl
from tme_datasets import list_registered_datasets, get_dataset_spec, Modality
from returns.maybe import Some

# 1. Retrieve the registered dataset collection
specs = list_registered_datasets()
print(f"Total registered datasets: {len(specs)}")  # e.g. 24

# 2. Dual indexing: by index or by dataset ID
first_spec = specs[0]
sade_spec = specs["GSE120575"]  # Direct ID indexing (raises KeyError if missing)
print(f"Dataset: {sade_spec.title} ({sade_spec.cancer_type})")

# 3. Attribute list accessor methods
all_ids = specs.ids()
all_titles = specs.titles()
all_modalities = specs.modalities()
organs_raw = specs.organs(unwrapped=True)  # Returns list[str | None]
sample_counts = specs.n_samples_or_cells(unwrapped=True)  # Returns list[int | None]

# 4. Safe lookup and functional filtering
match specs.get("GSE120575"):
    case Some(spec):
        print(f"Found: {spec.id} [{spec.platform}]")

sc_specs = specs.filter(modality=Modality.SINGLE_CELL)
print(f"Single-cell cohorts: {sc_specs.ids()}")

# 5. Convert to Polars DataFrame for downstream queries
df = specs.to_polars()
print(df.select(["id", "modality", "cancer_type", "platform", "has_response_labels"]))

# Run expressive Polars queries
melanoma_df = df.filter(pl.col("cancer_type") == "Melanoma")
print(f"Melanoma datasets: {melanoma_df['id'].to_list()}")
```

---

## 2. Loading Individual Datasets (Fast-Loading & H5AD Caching)

The `load_dataset` function serves as the central entry point for all single-cell references, bulk iAtlas validation cohorts, and direct paper datasets.

### Automatic Ingestion & Caching Workflow:
1. **Fast-Path (<0.5s)**: First checks if a preprocessed `.h5ad` file exists across configured paths (`data/preprocessed/{dataset_id}.h5ad`). If present, it loads directly from the binary H5AD format, completely bypassing slow gzipped text/tsv parsing.
2. **Auto-Download & Ingestion**: If no `.h5ad` file exists, it automatically downloads the vendor raw data (from NCBI GEO or cBioPortal) into `data/raw/{dataset_id}/` (if `auto_download=True`).
3. **Automatic H5AD Serialization**: Once parsed into an `AnnData` object, it immediately serializes the object to `data/preprocessed/{dataset_id}.h5ad`. All subsequent calls to `load_dataset` load from this H5AD file in milliseconds.
4. **Jupyter Notebook-Compatible Logging**: Logs stream immediately to `sys.stdout` with real-time download percentages, matrix parsing progress, and cache milestones without buffering delays or red stderr warnings.

```python
from tme_datasets import load_dataset
from returns.result import Success, Failure

# 1. Load single-cell dataset (loads from data/preprocessed/GSE120575.h5ad in <0.5s)
match load_dataset("GSE120575"):
    case Success(adata):
        print(f"Sade-Feldman AnnData: {adata.shape} (cells x genes)")
    case Failure(err):
        print(f"Failed to load: {err}")

# 2. Load bulk cohort (auto-downloads raw cBioPortal archive if missing, caches H5AD)
match load_dataset("Hugo-iAtlas"):
    case Success(adata):
        print(f"Hugo iAtlas: {adata.shape} (samples x genes)")
        print(f"Response labels present: {'response_binary' in adata.obs.columns}")

# 3. Force re-parsing and re-writing H5AD cache (e.g. after code update)
match load_dataset("GSE120575", force_recompute=True):
    case Success(adata):
        print("Recomputed and updated H5AD cache.")

# 4. Ensembl Normalization: automatically convert gene symbols to canonical Ensembl IDs
match load_dataset("GSE120575", normalize_ensembl=True, ensembl_release=111):
    case Success(adata):
        print(f"Ensembl normalized: {adata.shape}")
        print(f"Sample gene IDs: {list(adata.var_names[:3])}")
        print(f"Var attributes: {list(adata.var.columns)}")

# 5. Zero Hardcoding: Query dataset paths using the paths API
from tme_datasets.paths import (
    get_preprocessed_h5ad_path,
    get_raw_dataset_dir,
    get_ensembl_dir,
    find_dataset_h5ad,
)

h5ad_file = get_preprocessed_h5ad_path("GSE120575")
raw_dir = get_raw_dataset_dir("GSE120575")
ensembl_dir = get_ensembl_dir()
print(f"Canonical H5AD path:  {h5ad_file}")
print(f"Raw data directory:   {raw_dir}")
print(f"Ensembl cache folder: {ensembl_dir}")
```

---

## 3. Gene Set Collections & Signature Scoring

`tme_datasets.genesets` provides pre-registered TME signatures, custom GMT parsing, and signature scoring returning **Polars** DataFrames:

```python
from tme_datasets import (
    get_tme_major_lineage_collection,
    score_geneset_zscore,
    score_geneset_auc,
    compute_geneset_overlap,
    export_gmt,
    parse_gmt,
)

# 1. Load curated TME Major Lineage signatures
coll = get_tme_major_lineage_collection()
print(f"Collection: {coll.name}")
for name, gs in coll.gene_sets.items():
    print(f"  - {name}: {gs.genes[:5]}...")

# 2. Score signatures on an AnnData object
# Mean Z-Score
z_scores_df = score_geneset_zscore(adata, coll).unwrap()
print("Z-Score Signatures (Polars DataFrame):")
print(z_scores_df.head(3))

# Rank-based AUCell
auc_scores_df = score_geneset_auc(adata, coll, top_fraction=0.20).unwrap()
print("AUCell Rank Scores (Polars DataFrame):")
print(auc_scores_df.head(3))

# 3. Analyze overlap & redundancy between signatures
overlap_df = compute_geneset_overlap(coll, coll).unwrap()
print(overlap_df.filter(overlap_df["jaccard_similarity"] > 0).head())
```

---

## 4. Subsampling, Supersampling & Perturbations

Stochastic operations are implemented as pure functions returning new `AnnData` objects without in-place mutation.

### A. Subsampling & Supersampling
```python
from tme_datasets import subsample_cells, supersample_cells, SubsampleSpec
from returns.maybe import Some

# Subsample 50 cells balanced across cell types
sub_spec = SubsampleSpec(n_or_fraction=50, stratify_by=Some("cell_type"), balanced=True, seed=Some(42))
sub_adata = subsample_cells(adata, sub_spec).unwrap()
print(f"Subsampled shape: {sub_adata.shape}")

# Supersample (bootstrap with replacement) to 500 cells
super_adata = supersample_cells(adata, n_target=500, seed=Some(42)).unwrap()
print(f"Supersampled shape: {super_adata.shape}")
```

### B. Randomization & Normalization Methods Comparison

The package provides a comprehensive suite of count randomizations, parameter inference engines, and variance-stabilizing normalizations. The table below summarizes their mathematical formulations, output layers, zero-count behavior, and primary use cases:

| Method | Type | Mathematical Basis | Zero-Count Behavior | Output Layer / Var | Primary Use Case |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Simple Negative Binomial** | Randomization | $\lambda_{c,g} \sim \text{Gamma}(1/\alpha, \alpha X_{c,g})$, $Y \sim \text{Poisson}(\lambda)$ | $X_{c,g} = 0 \implies Y = 0$ (Zeros strictly preserved) | `.layers["randomized_nb"]` | Fast on-the-fly PyTorch data augmentation & noise sensitivity analysis |
| **Method of Moments (MoM)** | Inference + Resampling | $\hat{\mu}_g = \frac{\sum_i X_{i,g}}{\sum_i s_i}$, $\hat{\alpha}_g = \frac{\hat{\sigma}_g^2 - \hat{\mu}_g}{\hat{\mu}_g^2}$ | Technical dropouts sample $>0$; biological zeros remain $0$ | `.layers["randomized_nb"]`, `.varm["nb_means"]`, `.varm["nb_dispersions"]` | Fast analytical parameter recovery and dropout recovery across clusters |
| **Maximum Likelihood (MLE)** | Inference + Resampling | Profile likelihood Newton-Raphson on digamma $\psi(y + \frac{1}{\alpha})$ | Technical dropouts sample $>0$; biological zeros remain $0$ | `.layers["randomized_nb"]`, `.varm["nb_means"]`, `.varm["nb_dispersions"]` | Statistically optimal parameter estimation for moderate-sized cohorts |
| **Empirical Bayes (EB)** | Inference + Resampling | Parametric trend $\alpha(\mu) = a_0 + \frac{a_1}{\mu}$ with shrinkage | Technical dropouts sample $>0$; biological zeros remain $0$ | `.layers["randomized_nb"]`, `.varm["nb_means"]`, `.varm["nb_dispersions"]` | DESeq2/edgeR-style dispersion shrinkage for high-sparsity / low-count datasets |
| **Sanity Resampling** | Inference + Resampling | $n_{gc} \sim \text{Poisson}(N_c \alpha_g e^{\delta_{gc}})$, $\delta \sim \mathcal{N}(0, v_g)$ | Technical dropouts sample $>0$; biological zeros remain $0$ | `.layers["randomized_nb"]`, `.varm["nb_means"]`, `.varm["nb_dispersions"]` | First-principles sampling from inferred Log-Normal Poisson transcription states |
| **Sanity Normalization** | Normalization / Denoising | Laplace approximation of marginal likelihood $P(\mathbf{n}_g \mid v_g)$ | Zero counts shrunk toward dataset mean with wide error bars | `.layers["sanity_ltq"]`, `.layers["sanity_error"]`, `.var["sanity_variance"]` | Rigorous, parameter-free expression estimation and cell-to-cell distance calculation |
| **Analytic Pearson Residuals** | Normalization / HVG Selection | Closed-form offset NB: $r_{c,g} = \frac{n_{c,g} - \mu_{c,g}}{\sqrt{\mu_{c,g} + \mu_{c,g}^2 / \theta}}$ | Bounded negative residuals for zeros; no artificial zero-inflation | `.layers["pearson_residuals"]`, `.var["highly_variable"]` | Scalable variance stabilization and highly variable gene (HVG) selection |
| **Regularized GLM (sctransform)** | Normalization | $\log(\mu) \sim \beta_0 + \beta_1 \log_{10}(N)$, kernel smoothing across genes | Bounded negative residuals for zeros | `.layers["pearson_residuals"]`, `.var["sct_beta0"]`, `.var["sct_beta1"]` | Faithful port of Hafemeister & Satija (2019) with explicit sequencing-depth regularization |
| **Library Size CPM + $\log(1+x)$** | Normalization | $x'_{c,g} = \log\left(1 + 10^4 \cdot \frac{x_{c,g}}{N_c}\right)$ | Preserves exact zeros ($0 \to 0$) | `.X` (or custom layer) | Standard baseline preprocessing for downstream compatibility |
| **Capture Dropout Simulation** | Perturbation | Bernoulli mask: $P(\text{mask}_{c,g} = 0) = \text{rate}$ | Turns positive counts into exact zeros | `.X` | Simulating variable single-cell sequencing depth and capture inefficiencies |
| **Log-Normal Jitter** | Perturbation | Multiplicative noise: $X'_{c,g} = X_{c,g} \cdot e^{\mathcal{N}(0, \sigma^2)}$ | Preserves exact zeros | `.X` | Testing classifier and clustering stability against transcriptional noise |

---

### C. Detailed Usage of Each Method

#### 1. Simple Entry-Wise Negative Binomial Perturbation (Default)
Applies stochastic noise directly to each cell's observed expression level $X_{c,g}$ with fixed dispersion $\alpha$. If $X_{c,g} = 0$, it strictly yields $0$:
```python
from tme_datasets import randomize_negative_binomial, NegativeBinomialConfig
from returns.maybe import Some

# Default: Simple mode (zeros remain 0, unperturbed counts saved in .layers["raw_counts"])
nb_config = NegativeBinomialConfig(dispersion=0.20, seed=Some(123))
nb_adata = randomize_negative_binomial(adata, nb_config).unwrap()

print(f"Original counts preserved in: {list(nb_adata.layers.keys())}")
print(f"Mean count (raw): {adata.X.mean():.2f} | Mean count (NB): {nb_adata.X.mean():.2f}")
```

#### 2. Negative Binomial Parameter Inference & Resampling (MoM, MLE, Empirical Bayes)
Infers true underlying gene parameters $(\mu_g, \alpha_g)$ across cells (or stratified by `cluster_key="cell_type"`), accounting for cell library size depth $s_i$. Technical dropout zeros will sample counts $> 0$ with realistic probability, while true biological zeros remain 0:

```python
from tme_datasets import (
    randomize_negative_binomial,
    NegativeBinomialConfig,
    NBEstimationMethod,
    fit_nb_moments,
    fit_nb_mle,
    fit_nb_empirical_bayes,
)
from returns.maybe import Some

# Option A: Method of Moments (MoM) stratified by cell type
mom_config = NegativeBinomialConfig(
    estimation_method=Some(NBEstimationMethod.MOMENTS),
    cluster_key=Some("cell_type"),
    seed=Some(42),
)
mom_adata = randomize_negative_binomial(adata, mom_config).unwrap()
print(f"Inferred MoM means shape: {mom_adata.varm['nb_means'].shape}")

# Option B: Maximum Likelihood Estimation (MLE with Newton-Raphson on digamma)
mle_config = NegativeBinomialConfig(
    estimation_method=Some(NBEstimationMethod.MLE),
    seed=Some(42),
)
mle_adata = randomize_negative_binomial(adata, mle_config).unwrap()

# Option C: Empirical Bayes (Mean-dispersion trend fitting with shrinkage)
eb_config = NegativeBinomialConfig(
    estimation_method=Some(NBEstimationMethod.EMPIRICAL_BAYES),
    seed=Some(42),
)
eb_adata = randomize_negative_binomial(adata, eb_config).unwrap()

# Option D: Sanity Bayesian Log-Normal Poisson Resampling
sanity_nb_cfg = NegativeBinomialConfig(
    estimation_method=Some(NBEstimationMethod.SANITY),
    seed=Some(42),
)
sanity_resampled = randomize_negative_binomial(adata, sanity_nb_cfg).unwrap()
```

#### 3. Sanity Bayesian Denoising & Normalization (First-Principles LTQ Estimation)
Directly estimate Log-Transcription Quotients (LTQs) and analytical cell-by-gene error bars directly from raw counts (*Breda et al., Nature Biotechnology 2021*):

```python
from tme_datasets import run_sanity_normalization, SanityConfig

sanity_cfg = SanityConfig(v_min=0.001, v_max=20.0, n_bins=40)
sanity_adata = run_sanity_normalization(adata, sanity_cfg).unwrap()

print("Sanity Layers and Inferred Properties:")
print(f"  - Inferred LTQ shape:      {sanity_adata.layers['sanity_ltq'].shape}")
print(f"  - Posterior Error shape:    {sanity_adata.layers['sanity_error'].shape}")
print(f"  - Inferred Gene Variances:  {sanity_adata.var['sanity_variance'].head(3).to_dict()}")
```

#### 4. SCTransform & Analytic Pearson Residuals
Stabilize technical variance across sequencing depths using regularized Negative Binomial regression (*Hafemeister & Satija 2019*) or Analytic Pearson Residuals (*Lause et al. 2021*):

```python
from tme_datasets import normalize_sctransform, SCTransformConfig, SCTransformFlavor
from returns.maybe import Some

# Option A: Analytic Pearson Residuals (Fast, closed-form, native Scanpy)
analytic_cfg = SCTransformConfig(
    flavor=SCTransformFlavor.ANALYTIC,
    n_top_genes=Some(2000),  # Select top 2,000 HVGs via Pearson residual variance
    clip_residuals=True,     # Clip outliers to sqrt(N_cells)
)
sct_analytic_adata = normalize_sctransform(adata, analytic_cfg).unwrap()
print(f"Analytic Pearson Residuals: {sct_analytic_adata.layers['pearson_residuals'].shape}")
print(f"Top 2,000 HVGs selected in: .var['highly_variable']")

# Option B: Regularized GLM (Kernel-smoothed parameter regression)
glm_cfg = SCTransformConfig(
    flavor=SCTransformFlavor.REGULARIZED_GLM,
    clip_residuals=True,
)
sct_glm_adata = normalize_sctransform(adata, glm_cfg).unwrap()
print(f"Regularized GLM coefficients: .var['sct_beta0'], .var['sct_beta1']")
```

#### 5. Standard Library Size CPM + $\log(1+x)$ Normalization
```python
from tme_datasets import normalize_total_counts, log1p_transform

# Pure functional pipeline: Raw Counts -> CPM -> log(1 + CPM)
norm_adata = normalize_total_counts(adata, target_sum=1e4).unwrap()
log_adata = log1p_transform(norm_adata).unwrap()
print(f"CPM + log1p normalized max: {log_adata.X.max():.2f}")
```

### C. In-Silico Targeted Gene Knockout & Overexpression
Simulate targeted drug interventions or gene knockouts:

```python
from tme_datasets import in_silico_knockout, in_silico_overexpression

# Knockout checkpoint inhibitory receptors (PDCD1, CD274) with 100% efficiency
ko_adata = in_silico_knockout(adata, genes=("PDCD1", "CD274"), efficiency=1.0).unwrap()

# Overexpress IFNG by 3.0x fold-change
oe_adata = in_silico_overexpression(adata, genes=("IFNG",), fold_change=3.0).unwrap()
```

### D. Monadic Pipeline Composition
Chain multiple perturbations using `ComposeTransforms`:

```python
from tme_datasets import ComposeTransforms, simulate_dropout, add_expression_jitter

pipeline = ComposeTransforms([
    lambda a: subsample_cells(a, SubsampleSpec(n_or_fraction=0.8, seed=Some(1))),
    lambda a: simulate_dropout(a, rate=0.05, seed=2),
    lambda a: add_expression_jitter(a, sigma=0.05, seed=3),
])

perturbed_adata = pipeline(adata).unwrap()
```

---

## 5. In-Silico Pseudobulk Deconvolution Simulator

Generate synthetic bulk mixtures with exact known Dirichlet mixture proportions for benchmarking deconvolution algorithms (BayesPrism, InstaPrism):

```python
from tme_datasets import simulate_pseudobulk, PseudobulkConfig
from returns.maybe import Some

pb_config = PseudobulkConfig(
    n_samples=20,                 # 20 bulk RNA-seq mixtures
    cells_per_sample=1000,        # 1,000 single cells per mixture
    noise_dispersion=Some(0.05),  # Optional sequencing noise
    seed=Some(42),
)

bulk_adata, ground_truth_df = simulate_pseudobulk(
    adata,
    pb_config,
    cell_type_key="cell_type",
).unwrap()

print(f"Simulated Bulk Mixtures: {bulk_adata.shape} (Samples x Genes)")
print("Ground-Truth Proportions:")
print(ground_truth_df.head(3))
```

---

## 6. Multi-Dataset Harmonization & Integration Metrics

Query multiple datasets at once and automatically harmonize their feature and metadata spaces:

```python
from tme_datasets import query_datasets, HarmonizeConfig, HarmonizeMode, evaluate_integration_metrics

# Mode 1: Strict Gene Intersection
cfg_inter = HarmonizeConfig(mode=HarmonizeMode.INTERSECTION, min_shared_genes=500)
inter_adata = query_datasets(["Hugo-iAtlas", "Riaz-iAtlas"], config=cfg_inter).unwrap()
print(f"Intersection shape: {inter_adata.shape}")

# Mode 2: Union Zero-Filled (Pads unmeasured genes with zeros in sparse matrix)
cfg_union = HarmonizeConfig(mode=HarmonizeMode.UNION_ZERO_FILLED)
union_adata = query_datasets(["Hugo-iAtlas", "Riaz-iAtlas"], config=cfg_union).unwrap()
print(f"Union Zero-Filled shape: {union_adata.shape}")

# Evaluate Integration Quality & Mixing Metrics
metrics = evaluate_integration_metrics(
    inter_adata,
    batch_key="dataset_id",
    label_key="response_binary",
).unwrap()

print(f"Mean iLISI (Cohort Mixing, ideal -> 2.0):    {metrics.mean_ilisi}")
print(f"Mean cLISI (Cell-Type Purity, ideal -> 1.0): {metrics.mean_clisi}")
print(f"Silhouette Ratio (Bio / Batch):              {metrics.silhouette_ratio}")
print(f"kBET Acceptance Rate:                        {metrics.kbet_acceptance_rate}")
```

---

## 7. PyTorch Dataset & DataLoader Bridge

Train neural network models (autoencoders, scVI, classifier foundation models) with dynamic on-the-fly Negative Binomial augmentations:

```python
from tme_datasets import TmeTorchDataset, create_tme_dataloader
from returns.maybe import Some

# Initialize PyTorch Dataset
torch_ds = TmeTorchDataset(
    adata,
    label_keys=("response_binary", "cell_type"),
    on_the_fly_nb_dispersion=Some(0.15),  # Dynamic Gamma-Poisson sampling on each fetch
    device="cpu",
)

print(f"Total dataset items: {len(torch_ds)}")

# Create PyTorch DataLoader
loader = create_tme_dataloader(torch_ds, batch_size=32, shuffle=True)

for batch_x, batch_labels in loader:
    print(f"Batch X: {batch_x.shape} (Tensor float32)")
    print(f"Batch response labels: {batch_labels['response_binary'].shape}")
    break
```

---

## 8. Gene Identifier Unification & Ensembl Normalization

Heterogeneous datasets often arrive with varied identifier conventions: official HUGO gene symbols, previous symbols or synonyms, or Ensembl IDs with dot-version numbers (`ENSG00000133703.12`).

`tme_datasets` provides an industrial-strength, config-driven normalization engine that maps gene identifiers to canonical Ensembl gene IDs (defaulting to Ensembl Release 111, GRCh38), enriches `adata.var` with comprehensive genomic attributes, and resolves mapping conflicts deterministically.

### Key Capabilities:
- **PyEnsembl Auto-Installation**: Automatically downloads and indexes the required Ensembl release directly into `data/ensembl/` (never modifying user home directories).
- **Sub-Millisecond Persistent Parquet Caching**: Query results are stored in `data/ensembl/gene_mapping_cache_release_{release}.parquet` via **Polars**, reducing repeated normalization runs from minutes to under 50 milliseconds.
- **Deterministic Conflict Resolution**:
  - **1-to-many conflicts**: Prioritizes canonical reference chromosomes (`1-22`, `X`, `Y`, `MT`), `protein_coding` biotypes over pseudogenes, and selects the deterministic lowest alphanumeric ENSG ID. Unchosen IDs are stored in `adata.var["alternative_ensembl_ids"]`.
  - **0-to-many conflicts**: Resolves historical gene symbols, deprecated symbols, and synonyms via HGNC / `mygene` querying.
  - **Many-to-1 conflicts (duplicate Ensembl IDs)**: Expression values are aggregated using `sum`, `mean`, or `max`.
- **Rich `.var` Genomic Annotations**: Enriches every feature with:
  - `gene_id`: Canonical Ensembl gene ID (e.g. `ENSG00000133703`)
  - `gene_name`: Canonical gene symbol (e.g. `KRAS`)
  - `original_id`: Identifier originally present in the raw dataset
  - `contig`: Chromosome or scaffold (e.g. `12`)
  - `start`, `end`: Genomic coordinates
  - `strand`: `+` or `-`
  - `biotype`: e.g. `protein_coding`, `lncRNA`
  - `ensembl_release`, `species`
  - `mapping_status`: `exact_id`, `exact_symbol`, `alias_resolved`, `unmapped`
  - `alternative_ensembl_ids`: Comma-separated list of secondary Ensembl IDs

### Standalone Normalization Example:

```python
from tme_datasets import normalize_dataset_to_ensembl
from returns.result import Success, Failure

# 1. Normalize an AnnData object to Ensembl Release 111
match normalize_dataset_to_ensembl(adata, release=111, drop_unmapped=True, aggregation="sum"):
    case Success(norm_adata):
        print(f"Normalized shape: {norm_adata.shape}")
        print(f"Var index name:   {norm_adata.var_names.name}")  # 'gene_id'
        
        # Inspect genomic metadata added to var
        var_sample = norm_adata.var[["gene_name", "contig", "biotype", "mapping_status"]].head(5)
        print(var_sample)
        
        # Check unmapped genes recorded in uns
        if "unmapped_genes" in norm_adata.uns:
            print(f"Unmapped genes count: {norm_adata.uns['n_unmapped_genes']}")
    case Failure(err):
        print(f"Normalization failed: {err}")
```

### Low-Level Ensembl Normalization Usage:

```python
from tme_datasets.preprocessing.gene_normalization import (
    normalize_genes_to_ensembl,
    ensure_ensembl_release_installed,
)
from tme_datasets.paths import get_ensembl_dir

# 1. Auto-install and index Ensembl release into data/ensembl/
ensembl = ensure_ensembl_release_installed(release=111, species="human", ensembl_dir=get_ensembl_dir())
print(f"Ensembl release {ensembl.release} ready at: {get_ensembl_dir()}")

# 2. Directly normalize AnnData
norm_adata = normalize_genes_to_ensembl(
    adata,
    release=111,
    species="human",
    ensembl_dir=get_ensembl_dir(),
    drop_unmapped=True,
    aggregation="sum",
)
```

### Legacy Harmonization via `reconcile_genes`:

```python
from tme_datasets import reconcile_genes, GeneReconcileConfig, GeneIDType

# Quick symbol reconciliation (e.g. for cross-cohort intersection)
config = GeneReconcileConfig(
    target_type=GeneIDType.HUGO_SYMBOL,
    strip_version_suffix=True,
    handle_duplicates="sum",
)

reconciled_adata = reconcile_genes(adata, config).unwrap()
print(f"Reconciled feature count: {reconciled_adata.n_vars}")
```

---

## 9. Out-of-Core Storage & Cryptographic Verification

Efficiently manage large cohorts (Qian, Cheng, Tietscher) and verify data integrity:

```python
from tme_datasets import (
    compute_file_hash,
    verify_checksum,
    ChecksumSpec,
    load_backed,
    convert_to_zarr,
    StorageBackend,
)
from pathlib import Path

h5ad_path = Path("data/preprocessed/GSE120575.h5ad")

if h5ad_path.exists():
    # 1. Cryptographic SHA-256 Checksum
    sha256 = compute_file_hash(h5ad_path).unwrap()
    print(f"SHA-256: {sha256}")
    is_valid = verify_checksum(h5ad_path, ChecksumSpec(expected_hash=sha256)).unwrap()
    print(f"File verified: {is_valid}")

    # 2. Backed Mode (Read without RAM loading)
    backed_adata = load_backed(h5ad_path, backend=StorageBackend.BACKED_H5AD).unwrap()
    print(f"Backed AnnData: {backed_adata.shape} (is_backed={backed_adata.isbacked})")

    # 3. Convert to Chunked Zarr for high-performance streaming
    zarr_path = Path("data/preprocessed/GSE120575.zarr")
    convert_to_zarr(backed_adata.to_memory(), zarr_path).unwrap()
```

---

## 10. End-to-End Single-Cell Downstream Workflow (Case Study: GSE120575)

This section demonstrates a complete, production-grade downstream analysis pipeline for **Sade-Feldman et al. (Cell 2018, GSE120575)**: 16,291 CD45+ tumor-infiltrating immune cells from 48 metastatic melanoma patients treated with immune checkpoint blockade (ICB).

### Overview of Workflow Steps
1. **Ingestion & Caching**: Auto-download from NCBI GEO and cache as preprocessed H5AD to turn slow gzipped text parsing into an instantaneous load.
2. **Clinical Feature Engineering**: Standardize response labels (Responder vs Non-responder), patient ID, and treatment timing (Pre-baseline vs Post-treatment).
3. **Confounding Gene Filtering**: Purge mitochondrial (`MT-`), ribosomal (`RPS`/`RPL`), non-coding, and HLA artifacts.
4. **Variance Stabilization & HVG Selection**: Apply Analytic Pearson Residuals (`normalize_sctransform`) to compute exact residuals and identify the top 2,000 variable genes.
5. **Dimensionality Reduction & Graph Clustering**: PCA, kNN graph construction, UMAP embedding, and Leiden clustering.
6. **TME Lineage Annotation via Signature Scoring**: Score pre-registered TME lineage markers (`get_tme_major_lineage_collection()`) to assign CD8+ T, CD4+ T, B, NK, Monocyte/Macrophage, and Dendritic cell identities.
7. **Responder vs. Non-Responder Differential Abundance**: Calculate cell type composition differences across response groups using **Polars**.
8. **In-Silico Checkpoint Knockout & Pseudobulk Simulation**: Perform in-silico target ablation (*PDCD1*, *CTLA4*) and generate ground-truth pseudobulk mixtures for deconvolution benchmarks.

---

### Step 1: Ingestion & Fast Local Caching

Loading directly from raw GEO text files can take several minutes to decompress. Caching the parsed `AnnData` to H5AD once allows all subsequent sessions to load in milliseconds:

```python
from pathlib import Path
from tme_datasets import load_dataset
from returns.result import Success, Failure
import anndata as ad

cache_file = Path("data/preprocessed/GSE120575.h5ad")

if cache_file.exists():
    print(f"Loading cached AnnData from {cache_file}...")
    adata = ad.read_h5ad(cache_file)
else:
    print("Loading GSE120575 via tme_datasets (auto-downloads from GEO if missing)...")
    match load_dataset("GSE120575"):
        case Success(raw_adata):
            adata = raw_adata
            cache_file.parent.mkdir(parents=True, exist_ok=True)
            adata.write_h5ad(cache_file)
            print(f"Saved cached dataset to {cache_file}")
        case Failure(err):
            raise RuntimeError(f"Could not load GSE120575: {err}")

print(f"Initial dimensions: {adata.shape} (16,291 cells x 55,737 genes)")
```

---

### Step 2: Clinical Feature Engineering & Standardization

The raw metadata contains columns such as `characteristics: response` and `characteristics: patinet ID (Pre=baseline; Post= on treatment)`. We standardize these into clean observation columns:

```python
import polars as pl
import pandas as pd

# 1. Binarize response (Responder -> 1, Non-responder -> 0)
if "characteristics: response" in adata.obs.columns:
    resp_map = {"Responder": 1, "Non-responder": 0}
    adata.obs["response_binary"] = adata.obs["characteristics: response"].map(resp_map)
    print(f"Response breakdown:\n{adata.obs['response_binary'].value_counts(dropna=False)}")

# 2. Extract patient identifier and timepoint (Pre vs Post)
if "characteristics: patinet ID (Pre=baseline; Post= on treatment)" in adata.obs.columns:
    pat_series = adata.obs["characteristics: patinet ID (Pre=baseline; Post= on treatment)"].astype(str)
    
    # Extract timepoint: Pre (baseline) vs Post (on-treatment)
    adata.obs["timepoint"] = pat_series.apply(
        lambda s: "Pre" if "Pre" in s else ("Post" if "Post" in s else "Unknown")
    )
    # Extract clean patient ID: e.g. P1, P2, P33
    adata.obs["patient_id"] = pat_series.str.extract(r"(P\d+)")

# 3. Therapy regimen (anti-PD1, anti-CTLA4, anti-CTLA4+PD1)
if "characteristics: therapy" in adata.obs.columns:
    adata.obs["therapy"] = adata.obs["characteristics: therapy"]

print(f"Timepoints: {dict(adata.obs['timepoint'].value_counts())}")
print(f"Therapies:  {dict(adata.obs['therapy'].value_counts())}")
```

---

### Step 3: Confounding Gene Filtering

Single-cell immune profiling is frequently contaminated by non-informative mitochondrial, ribosomal, and immunoglobulin/HLA artifacts. We filter these out using `filter_confounding_genes`:

```python
from tme_datasets import filter_confounding_genes

# Remove MT-, RPS-, RPL-, HLA-, and non-coding RNA confounds
clean_adata = filter_confounding_genes(adata).unwrap()
print(f"Post-filtering dimensions: {clean_adata.shape} (removed {adata.n_vars - clean_adata.n_vars} confounding genes)")
```

---

### Step 4: Variance Stabilization (Analytic Pearson Residuals)

We apply **Analytic Pearson Residuals** (*Lause et al., Genome Biology 2021*) to stabilize variance, rank the top 2,000 highly variable genes, and bound outlier residuals:

```python
from tme_datasets import normalize_sctransform, SCTransformConfig, SCTransformFlavor
from returns.maybe import Some

sct_cfg = SCTransformConfig(
    flavor=SCTransformFlavor.ANALYTIC,
    n_top_genes=Some(2000),
    clip_residuals=True,
    use_layer_as_x=False,  # Preserves raw counts in .X and stores residuals in .layers['pearson_residuals']
)

norm_adata = normalize_sctransform(clean_adata, sct_cfg).unwrap()
print(f"Pearson residuals layer shape: {norm_adata.layers['pearson_residuals'].shape}")
print(f"Highly variable genes selected: {norm_adata.var['highly_variable'].sum()}")
```

> [!TIP]
> Alternatively, for true parameter-free Bayesian expression states and posterior standard error bars (without library size heuristics), use `run_sanity_normalization(clean_adata)` (*Breda et al., Nature Biotechnology 2021*).

---

### Step 5: Dimensionality Reduction, Neighborhood Graph & Clustering

Using Scanpy on the variance-stabilized Pearson residuals:

```python
import scanpy as sc

# Work on Pearson residuals for PCA and graph construction
processed_adata = norm_adata.copy()
processed_adata.X = processed_adata.layers["pearson_residuals"].copy()

# 1. Principal Component Analysis on HVGs
sc.tl.pca(processed_adata, n_comps=30, use_highly_variable=True)

# 2. k-Nearest Neighbors Graph
sc.pp.neighbors(processed_adata, n_neighbors=15, n_pcs=30)

# 3. UMAP Embedding
sc.tl.umap(processed_adata)

# 4. Leiden Community Detection
sc.tl.leiden(processed_adata, resolution=0.5, key_added="leiden_cluster")
print(f"Identified {processed_adata.obs['leiden_cluster'].nunique()} Leiden clusters.")
```

---

### Step 6: TME Lineage Annotation via Signature Scoring

We score each cell against pre-registered tumor microenvironment lineage marker collections (`get_tme_major_lineage_collection()`), which returns a **Polars** DataFrame with standardized Z-scores or AUC ranks:

```python
from tme_datasets import get_tme_major_lineage_collection, score_geneset_zscore
import numpy as np

# Retrieve curated TME major lineage gene sets:
# CD8_T_cell, CD4_T_cell, B_cell, NK_cell, Monocyte_Macrophage, Dendritic_cell, etc.
lineage_coll = get_tme_major_lineage_collection()
print(f"Scoring {len(lineage_coll.gene_sets)} lineage signatures across cells...")

# Score signatures using standardized Z-score
scores_df = score_geneset_zscore(norm_adata, lineage_coll).unwrap()

# Join scores back into AnnData observation metadata
for col in scores_df.columns:
    if col not in ("cell_id", "sample_id"):
        processed_adata.obs[f"sig_{col}"] = scores_df[col].to_numpy()

# Annotate clusters based on highest mean lineage score
cluster_lineage_map = {}
for cluster_id in processed_adata.obs["leiden_cluster"].unique():
    cluster_cells = processed_adata.obs["leiden_cluster"] == cluster_id
    best_lineage = None
    best_score = -np.inf
    for gs in lineage_coll.gene_sets.values():
        mean_score = processed_adata.obs.loc[cluster_cells, f"sig_{gs.name}"].mean()
        if mean_score > best_score:
            best_score = mean_score
            best_lineage = gs.name
    cluster_lineage_map[cluster_id] = best_lineage

processed_adata.obs["cell_type"] = processed_adata.obs["leiden_cluster"].map(cluster_lineage_map)
print("Annotated Cell-Type Distribution:")
print(processed_adata.obs["cell_type"].value_counts())
```

---

### Step 7: Responder vs. Non-Responder Differential Composition Analysis

Using **Polars**, compute the relative cell-type fractions across response groups to examine immune microenvironment remodeling:

```python
import polars as pl

# Convert metadata to Polars DataFrame
obs_df = pl.from_pandas(processed_adata.obs.reset_index())

# Group by clinical response and cell type
composition = (
    obs_df.filter(pl.col("response_binary").is_not_null())
    .group_by(["response_binary", "cell_type"])
    .len()
    .with_columns(
        (pl.col("len") / pl.col("len").sum().over("response_binary")).alias("fraction")
    )
    .sort(["cell_type", "response_binary"])
)

print("\n--- Cell-Type Proportions: Responders (1) vs Non-Responders (0) ---")
print(composition)
```

---

### Step 8: In-Silico Checkpoint Perturbation & Pseudobulk Simulation

#### A. In-Silico Knockout of Immune Checkpoint Targets
Simulate the functional loss or blockade of inhibitory receptors (*PDCD1*, *CTLA4*, *HAVCR2*, *LAG3*):

```python
from tme_datasets import in_silico_knockout, PerturbationConfig

ko_cfg = PerturbationConfig(
    target_genes=("PDCD1", "CTLA4", "HAVCR2"),
)

ko_adata = in_silico_knockout(processed_adata, ko_cfg).unwrap()
print("Knockout verified: expression of target genes set strictly to zero.")
```

#### B. Generate Synthetic Bulk Mixtures with Exact Ground Truth
Simulate 50 patient pseudobulk RNA-seq mixtures (1,000 cells per sample) with Dirichlet-sampled cell proportions to benchmark deconvolution algorithms (*BayesPrism*, *InstaPrism*):

```python
from tme_datasets import simulate_pseudobulk, PseudobulkConfig
from returns.maybe import Some

bulk_cfg = PseudobulkConfig(
    n_samples=50,
    cells_per_sample=1000,
    noise_dispersion=Some(0.1),  # Add biological overdispersion noise
    seed=Some(42),
)

# Simulate mixtures with exact known proportions
sim_bulk_adata, truth_proportions_df = simulate_pseudobulk(
    processed_adata,
    bulk_cfg,
    cell_type_key="cell_type",
).unwrap()

print(f"Simulated bulk matrix: {sim_bulk_adata.shape} (50 bulk samples x genes)")
print("Ground truth cell fractions (first 3 samples):")
print(truth_proportions_df.head(3))
```

---

## Summary

`tme_datasets` provides a unified, mathematically pure, and test-covered functional core for all dataset querying, harmonization, and simulation needs across the project.
