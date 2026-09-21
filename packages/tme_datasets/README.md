# `tme_datasets`

A unified data access, downloading, preprocessing, harmonization, and perturbation framework for all single-cell reference atlases, spatial transcriptomics, bulk clinical validation cohorts, and gene set collections across tumor microenvironment analyses.

This package consolidates and replaces:
- `single_cell_datasets`
- `single_cell_immuno_datasets`
- `genentech_datasets`
- `ici_datasets`

---

## Key Features

1. **Declarative Query API**:
   - Query individual single-cell or bulk cohorts: `load_dataset("GSE120575")`, `load_dataset("Hugo-iAtlas")`, `load_dataset("Auslander")`.
   - Multi-dataset harmonization: `query_datasets(["Hugo-iAtlas", "Riaz-iAtlas"], config=HarmonizeConfig(mode=HarmonizeMode.INTERSECTION))`.
   - Supports both `INTERSECTION` and `UNION_ZERO_FILLED` (sparse padded) feature alignment.

2. **Gene Set Datasets & Signature Scoring (`tme_datasets.genesets`)**:
   - Curated collections: Bagaev MFP (29 signatures), ImmunoCompass, MSigDB Hallmarks/C7, TME lineage markers, Ayers et al. 18-gene GEP.
   - Built-in scoring: `score_geneset_zscore`, `score_geneset_auc` (AUCell), and `compute_geneset_overlap`.

3. **Stochastic Perturbations & Sampling (`tme_datasets.transforms`)**:
   - `subsample_cells`: Uniform, stratified, or class-balanced downsampling.
   - `supersample_cells`: Bootstrap oversampling with replacement.
   - `randomize_negative_binomial`: Gamma-Poisson mixture model simulating realistic single-cell biological overdispersion.
   - `simulate_dropout` and `add_expression_jitter`.
   - `in_silico_knockout` and `in_silico_overexpression`.
   - `ComposeTransforms`: PyTorch-style pipeline chaining.

4. **In-Silico Pseudobulk Simulator (`tme_datasets.simulation`)**:
   - `simulate_pseudobulk()` generates synthetic bulk mixtures with exact known Dirichlet mixture proportions.

5. **PyTorch Dataset & DataLoader Bridge (`tme_datasets.torch`)**:
   - `TmeTorchDataset` delivers batches with on-the-fly Negative Binomial augmentations directly for deep learning training.
   - `create_tme_dataloader`.

6. **Integration Quality & Mixing Metrics (`tme_datasets.harmonization`)**:
   - Quantitative batch assessment: `iLISI`, `cLISI`, `batch_silhouette`, and `kBET`.

7. **Gene Identifier Reconciliation (`tme_datasets.genes`)**:
   - Automated conversion and collision resolution across Ensembl IDs, HGNC Symbols, and Entrez IDs.

8. **Out-of-Core Processing (`tme_datasets.storage`)**:
   - Backed AnnData (`backed="r"`) and chunked Zarr conversion for multi-hundred-thousand cell cohorts.

9. **Cryptographic Verification (`tme_datasets.download`)**:
   - SHA-256 chunked hashing for data integrity.

---

## Quickstart

```python
from tme_datasets import (
    load_dataset,
    query_datasets,
    HarmonizeConfig,
    HarmonizeMode,
    subsample_cells,
    SubsampleSpec,
    randomize_negative_binomial,
    NegativeBinomialConfig,
    load_geneset_collection,
    score_geneset_zscore,
)
from returns.result import Success, Failure
from returns.maybe import Some

# 1. Load an individual dataset
match load_dataset("GSE120575"):
    case Success(adata):
        print(f"Loaded Sade-Feldman: {adata.shape}")
    case Failure(err):
        print(f"Error: {err}")

# 2. Query and harmonize multiple bulk cohorts
match query_datasets(["Hugo-iAtlas", "Riaz-iAtlas"], config=HarmonizeConfig(mode=HarmonizeMode.INTERSECTION)):
    case Success(combined):
        print(f"Harmonized cohort: {combined.shape}")

# 3. Apply Negative Binomial perturbation
config = NegativeBinomialConfig(dispersion=0.15, seed=Some(42))
match randomize_negative_binomial(adata, config):
    case Success(perturbed):
        print("Perturbation complete")

# 4. Score TME gene signatures
match load_geneset_collection("tme_major_lineages"):
    case Success(coll):
        scores_df = score_geneset_zscore(adata, coll).unwrap()
        print(scores_df.head())
```
