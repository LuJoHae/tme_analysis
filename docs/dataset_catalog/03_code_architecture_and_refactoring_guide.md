# Code Architecture & Refactoring Guide for AI Agents

This guide provides concrete, actionable technical instructions for any AI agent or engineer refactoring, extending, or maintaining data pipelines in this codebase.

---

## 1. End-to-End Pipeline Execution Flow

```
scripts/01_download_datasets.py
    |
    v
scripts/sade_feldman_deconv_validation/
    |
    +--> 01_build_reference.py (Sade-Feldman standalone reference phi)
    |
    +--> 01b_build_integrated_reference.py (Harmony integration with 10x atlas)
    |
    +--> 01c_build_dataset_reference.py (Single dataset references: Jerby, Maynard, Ma, Yost)
    |
    +--> 01d_build_combined_references.py (Combinations of n datasets & cell subsamplings)
    |
    v
02_deconvolute_iatlas.py (BayesPrism / InstaPrism NNLS deconvolution of 9 bulk iAtlas cohorts)
    |
    v
03_logistic_regression.py (Univariate & multivariate logistic regression vs RECIST response)
    |
    v
04_analyze_milopy.py (Orthogonal single-cell Milo differential abundance on Sade-Feldman)
    |
    v
05_compare_concordance.py (Compiles multi_resolution_benchmark_summary.parquet)
    |
    v
06_plot_figures.py (Altair SVG/PNG vector generation: Fig 1-4 and Step 1-6 figures)
```

---

## 2. File-by-File Code Map: Data Ingestion & Processing

### A. Raw Dataset Downloaders & Loaders

#### 1. [`scripts/01_download_datasets.py`](file:///Users/halu/Code/tme_analysis/scripts/01_download_datasets.py)
* **Function**: `fetch_signature_dataset(name: str, base_dir: Path)`
* **Lines**: 97–155
* **Behavior**: Instantiates classes from `singlecellrnasignature.adata` (e.g., `PelkaSpatiallyOrganizedMulticellular2021Adata`), derives them via `datalair`, and writes clean `.h5ad` files to `base_dir / name / {name}.h5ad`.
* **Refactoring Note**: If downloading from GEO fails or times out, download cache in `~/.cache/datalair` and `/tmp` is wiped automatically.

#### 2. [`packages/single_cell_datasets/src/single_cell_datasets/_single_cell_datasets.py`](file:///Users/halu/Code/tme_analysis/packages/single_cell_datasets/src/single_cell_datasets/_single_cell_datasets.py)
* **Classes**: `SingleCellDataProcessStep01` through `Step07`
* **Lines**: 77–385
* **Behavior**: Defines the multi-step aggregation of 17 single-cell datasets into a unified pan-cancer atlas:
  - `Step01`: Ingests and standardizes each study.
  - `Step02`: Harmonizes metadata schemas (`patient`, `organ`, `cancer_type`, `is_tumor`).
  - `Step07`: Combines all 10x matrices and exports `adata.h5ad`.
* **Important Caveat for Agents**: Line 78 hardcodes `storage_path = Path("/storage/halu").resolve()`. If running on a local workstation without `/storage/halu`, use local paths or `scratch/lair`.

---

### B. Single-Cell Reference Generation

#### 1. [`scripts/sade_feldman_deconv_validation/01_build_reference.py`](file:///Users/halu/Code/tme_analysis/scripts/sade_feldman_deconv_validation/01_build_reference.py)
* **Input**: `scratch/GSE120575/gse120575_tpm.parquet` (16,288 cells).
* **Gene Filtering** (Lines 75–95): Excludes ribosomal (`RP[SL]`), mitochondrial (`MT-`), immunoglobulins (`IGH|IGK|IGL`), pseudogenes/lncRNAs (`RP11|AC00|RNU|RNA5S|LINC|CTD-`), and melanoma markers (`MLANA|PMEL|TYR|DCT|MITF|S100B|MAGEA`).
* **Output**:
  - `reference_phi.parquet` ($K \times G$ linear simplex matrix).
  - `reference_marker_genes.parquet` (Top 35 marker genes per cluster).

#### 2. [`scripts/sade_feldman_deconv_validation/01c_build_dataset_reference.py`](file:///Users/halu/Code/tme_analysis/scripts/sade_feldman_deconv_validation/01c_build_dataset_reference.py)
* **Standalone Loaders**:
  - `load_jerby_arnon(raw_dir)` (Lines 140–192): GSE115978 (7,186 cells).
  - `load_maynard(repo_root)` (Lines 194–224): Maynard NSCLC (3,000 cells).
  - `load_ma_liver(raw_dir)` (Lines 226–280): GSE125449 (5,115 cells).
  - `load_yost(raw_dir)` (Lines 282–340): GSE123813 (3,500 cells).
* **Multi-Resolution Clustering**: Runs Leiden clustering across 8 resolutions:
  `DEFAULT_RESOLUTIONS = (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0)`.
* **Condition Number**: Computes 2-norm condition number $\kappa(\Phi) = \frac{\sigma_{\max}(\Phi)}{\sigma_{\min}(\Phi)}$ to quantify collinearity.

#### 3. [`scripts/sade_feldman_deconv_validation/01d_build_combined_references.py`](file:///Users/halu/Code/tme_analysis/scripts/sade_feldman_deconv_validation/01d_build_combined_references.py)
* **Purpose**: Constructs combinations of $n$ datasets and evaluates cell subsampling ladders.
* **Component Loaders** (Lines 75–104): Caches datasets in memory (`dataset_cache`) so matrices are loaded only once.
* **Harmony Integration** (Lines 145–215):
  - Subsets shared gene intersection across component datasets (requires $\ge 500$ shared genes).
  - Corrects dataset-specific batch effects using `harmonypy.run_harmony`.
  - Fixes `harmonypy` PCA orientation: checks whether `ho.Z_corr` is shaped `(n_pcs, n_cells)` or `(n_cells, n_pcs)` and transposes accordingly.
* **Cell Subsampling Engine** (Lines 375–412):
  - For each combination, subsamples a fraction $\alpha \in [0.10, 1.00]$ of cells preserving cohort balance.
  - Builds reference signatures at resolution `0.5` and saves metrics to `reference_resolution_metrics_{combo_id}.parquet`.

---

### C. Bulk Deconvolution Engine

#### [`scripts/sade_feldman_deconv_validation/02_deconvolute_iatlas.py`](file:///Users/halu/Code/tme_analysis/scripts/sade_feldman_deconv_validation/02_deconvolute_iatlas.py)
* **Input References**: Scans `reference_phi_*.parquet`.
* **Input Bulk**: Scans `scratch/lair/CBioPortalDataset-*` (or `/storage/halu/lair`).
* **Caching Mechanism** (Lines 340–355): Checks if `deconv_fractions_{tag}.parquet` already exists and exceeds 1,000 bytes. Skips computation if cached.
* **Algorithm**:
  - `instaprism_deconvolve(X_bulk, Phi_ref)`: Linear Non-Negative Least Squares (NNLS) with simplex constraint $\sum_k \theta_k = 1$.
  - Intersects genes between bulk expression and reference $\Phi$.
  - Concatenates inferred fractions across all 9 cohorts ($N=1,097$ samples) into a single parquet.

---

### D. Response Prediction & Benchmark Aggregation

#### 1. [`scripts/sade_feldman_deconv_validation/03_logistic_regression.py`](file:///Users/halu/Code/tme_analysis/scripts/sade_feldman_deconv_validation/03_logistic_regression.py)
* **Input**: `deconv_fractions_{tag}.parquet`.
* **Clinical Association** (Lines 110–235):
  - Fits univariate logistic regression for each cell state: $\beta$, standard error, p-value, odds ratio, FDR (Benjamini-Hochberg), and univariate ROC-AUC.
  - Fits multivariate logistic regression (L2 regularization) on all cell states simultaneously to compute `multivariate_auc`.
  - Stratifies across `Melanoma` ($N=347$), `Pan-Cancer` ($N=1,097$), and individual cohorts.
* **Output**: `logistic_results_{tag}.parquet` and combined `logistic_regression_results_all_resolutions.parquet`.

#### 2. [`scripts/sade_feldman_deconv_validation/05_compare_concordance.py`](file:///Users/halu/Code/tme_analysis/scripts/sade_feldman_deconv_validation/05_compare_concordance.py)
* **Input**: All `logistic_results_*.parquet` and `reference_resolution_metrics_*.parquet`.
* **Benchmark Compilation** (Lines 360–540):
  - Aggregates metrics across all reference configurations.
  - Stores `reference_type`, `reference_category`, `resolution`, `n_clusters`, `n_signature_genes`, `condition_number`, `melanoma_multivariate_auc`, `pancancer_multivariate_auc`, `mean_cohort_multivariate_auc`, `n_datasets`, `n_cells`, `subsample_fraction`.
* **Output**: `multi_resolution_benchmark_summary.parquet`.

#### 3. [`scripts/sade_feldman_deconv_validation/06_plot_figures.py`](file:///Users/halu/Code/tme_analysis/scripts/sade_feldman_deconv_validation/06_plot_figures.py)
* **Figure 1 Generator**: `plot_fig1_predictive_capacity(data_dir, results_dir)` (Lines 1227–1615).
* **6-Row 18-Panel Architecture**:
  - Row 1 (A–C): Resolution vs. AUC (Melanoma, Pan-Cancer, Mean Cohort).
  - Row 2 (D–F): Cell States (Clusters) vs. AUC.
  - Row 3 (G–I): Number of Datasets ($n \in \{1, 2, 3, 4, 16\}$) vs. AUC.
  - Row 4 (J–L): Single Cells ($N_{\text{cells}} \in [1.5k, 41.3k]$) vs. AUC with logarithmic saturation curves.
  - Row 5 (M–O): Cell Subsampling Fraction ($\alpha \in [10\%, 100\%]$) titration curves.
  - Row 6 (P–R): Strategy Boxplots and Marginal Gain per Resolution Step.
* **Output**: `results/sade_feldman_deconv_validation/fig1_predictive_capacity_benchmark.svg` and `.png`.

---

## 3. How to Refactor for Full-Depth Datasets (Removing the 2,083 Cap)

To replace the 2,083 subsampled cohorts with full datasets, another AI agent should follow these exact steps:

### Step 1: Add Sade-Feldman to `01d_build_combined_references.py`
In `scripts/sade_feldman_deconv_validation/01d_build_combined_references.py`:
1. Import `load_sade_feldman`:
   ```python
   from scripts.sade_feldman_deconv_validation.01_build_reference import load_sade_feldman_raw
   ```
2. In `load_dataset_cached`, add the `"sf"` case:
   ```python
   case "sf":
       res = load_sade_feldman(raw_dir / "GSE120575")
   ```
3. Update `CRITERIA_COMBOS` and `RANDOM_COMBOS` to include quintuplets ($n=5$):
   ```python
   "all_5_combined": ("All-5-Datasets-Combined", "Criteria-Combined", ("sf", "jerby", "maynard", "ma", "yost")),
   ```
   Total cell depth will be $16{,}288 + 7{,}186 + 5{,}115 + 3{,}500 + 3{,}000 = \mathbf{35{,}089}\text{ authentic cells}$.

### Step 2: Ingesting Additional Full Single-Cell Datasets
To incorporate `PelkaSpatiallyOrganizedMulticellular2021` (~65k cells), `Azizi` (~45k cells), or `Qian` (~200k cells):
1. Download via `scripts/01_download_datasets.py`:
   ```bash
   .venv/bin/python scripts/01_download_datasets.py --out-dir data/raw
   ```
2. Add corresponding loaders in `01d_build_combined_references.py`:
   - Load `.h5ad` via `anndata.read_h5ad`.
   - Filter confounding gene families using `filter_confounding_genes`.
   - Compute linear counts/TPM matrix.
3. Add dataset IDs to `target_combos` in `01d_build_combined_references.py`.

### Step 3: Run the Cascade
Execute the pipeline sequentially:
```bash
# 1. Build combined reference signatures
.venv/bin/python scripts/sade_feldman_deconv_validation/01d_build_combined_references.py \
  --out-dir "output/output/sade_feldman_deconv_validation" \
  --subsample-combos "all_5_combined,random_triplet1"

# 2. Deconvolute bulk iAtlas cohorts
.venv/bin/python scripts/sade_feldman_deconv_validation/02_deconvolute_iatlas.py \
  --out-dir "output/output/sade_feldman_deconv_validation" \
  --benchmark-mode

# 3. Fit logistic regression models
.venv/bin/python scripts/sade_feldman_deconv_validation/03_logistic_regression.py \
  --out-dir "output/output/sade_feldman_deconv_validation" \
  --benchmark-mode

# 4. Aggregate benchmark summary
.venv/bin/python scripts/sade_feldman_deconv_validation/05_compare_concordance.py \
  --out-dir "output/output/sade_feldman_deconv_validation"

# 5. Render updated vector figures
.venv/bin/python scripts/sade_feldman_deconv_validation/06_plot_figures.py \
  --data-dir "output/output/sade_feldman_deconv_validation" \
  --results-dir "results/sade_feldman_deconv_validation"
```
