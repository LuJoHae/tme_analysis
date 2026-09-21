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

---

## 1. Dataset Registry & Metadata Inspection

All single-cell references, bulk iAtlas cohorts, and direct publication datasets are declaratively indexed in an immutable registry.

```python
from tme_datasets import list_registered_datasets, get_dataset_spec
from returns.maybe import Some

# List all available datasets across modalities
specs = list_registered_datasets()
print(f"Total registered datasets: {len(specs)}")

# Inspect a specific dataset's metadata specification
match get_dataset_spec("GSE120575"):
    case Some(spec):
        print(f"ID:       {spec.id}")
        print(f"Title:    {spec.title}")
        print(f"Modality: {spec.modality.value}")
        print(f"Cancer:   {spec.cancer_type}")
        print(f"Platform: {spec.platform}")
        print(f"Response: {spec.has_response_labels}")
```

---

## 2. Loading Individual Datasets

The `load_dataset` function serves as a single entry point for all datasets, automatically resolving local paths and parsing clinical metadata:

```python
from tme_datasets import load_dataset
from returns.result import Success, Failure

# Load a single-cell dataset
match load_dataset("GSE120575"):
    case Success(adata):
        print(f"Sade-Feldman AnnData: {adata.shape} (cells x genes)")
    case Failure(err):
        print(f"Failed to load: {err}")

# Load a bulk immunotherapy cohort (cBioPortal / iAtlas track)
match load_dataset("Hugo-iAtlas"):
    case Success(adata):
        print(f"Hugo iAtlas: {adata.shape} (samples x genes)")
        print(f"Response labels present: {'response_binary' in adata.obs.columns}")

# Load a direct paper dataset (Publication track)
match load_dataset("Auslander"):
    case Success(adata):
        print(f"Auslander direct paper: {adata.shape}")
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

## 8. Gene Identifier Unification

Reconcile heterogeneous gene identifiers (Ensembl version suffixes, HGNC symbols, Entrez):

```python
from tme_datasets import reconcile_genes, GeneReconcileConfig, GeneIDType

# Automatically strips ENSG...15 suffixes and maps to HUGO symbols
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

## Summary

`tme_datasets` provides a unified, mathematically pure, and test-covered functional core for all dataset querying, harmonization, and simulation needs across the project.
