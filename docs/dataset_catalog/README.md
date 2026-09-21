# Dataset Catalog & Reference Architecture Guide

This directory provides an exhaustive, authoritative catalog of all single-cell transcriptomic reference datasets and bulk clinical validation cohorts used across this repository (`tme_analysis`). It is structured specifically so that subsequent AI agents and human researchers can understand the data provenance, navigate the code, and refactor data ingestion pipelines with zero ambiguity.

---

## 1. High-Level Dataset Classification

The project operates at the interface of **single-cell reference atlases** and **bulk clinical immunotherapy cohorts**:

```
+----------------------------------------------------------------------------------------------------+
|                                     SINGLE-CELL REFERENCE POOL                                     |
|  Purpose: Define immune cell states & calculate linear cell-state expression signatures (Phi).     |
|  Response Labels: NOT REQUIRED (except for orthogonal Milo single-cell DA testing in Sade-Feldman) |
+----------------------------------------------------------------------------------------------------+
              |
              | 1. Gene Filtering (remove ribosomal, MT, BCR/TCR, tumor antigens)
              | 2. Harmony Batch Integration across platforms & cohorts
              | 3. Multi-Resolution Leiden Clustering (res in [0.25, 2.00])
              | 4. Linear Simplex Cell-State Signature (Phi in R^{K x G}, sum_g Phi_{c,g} = 1)
              v
+----------------------------------------------------------------------------------------------------+
|                                       BULK DECONVOLUTION ENGINE                                    |
|  BayesPrism / InstaPrism non-negative least squares (NNLS) deconvolution on bulk RNA-seq cohorts   |
+----------------------------------------------------------------------------------------------------+
              ^
              | 5. Bulk Mixture Expression (X in R^{G x S})
              |
+----------------------------------------------------------------------------------------------------+
|                                   BULK CLINICAL VALIDATION COHORTS                                 |
|  Purpose: Deconvolute bulk mixtures and test if cell fractions predict clinical response.          |
|  Response Labels: MANDATORY (RECIST Responder vs. Non-Responder; OS / PFS survival endpoints)      |
+----------------------------------------------------------------------------------------------------+
              |
              v
+----------------------------------------------------------------------------------------------------+
|                                  BENCHMARK EVALUATION & DISCOVERY                                  |
|  - Logistic Regression (Univariate & Multivariate ROC-AUC in Melanoma & Pan-Cancer strata)         |
|  - Cross-Modality Concordance (Bulk Regression Betas vs. Single-Cell Milo log2FC)                  |
|  - Reference Scaling Benchmarks (Predictive AUC vs. Resolution, Clusters, Datasets, Cells)         |
+----------------------------------------------------------------------------------------------------+
```

---

## 2. Catalog Directory Structure

| File | Scope & Contents |
| :--- | :--- |
| [`01_single_cell_reference_datasets.md`](./01_single_cell_reference_datasets.md) | Exhaustive documentation of all single-cell datasets (Sade-Feldman, Jerby-Arnon, Ma, Yost, Maynard, GSE179994, and all 17 pan-cancer atlas cohorts), detailing biological origin, cell counts, file locations, preprocessing, and the 2,083-cell subsampling origin. |
| [`02_bulk_validation_cohorts.md`](./02_bulk_validation_cohorts.md) | Comprehensive catalog of the 9 cBioPortal / iAtlas bulk cohorts (1,097 samples) and historical bulk cohorts, including sample sizes, cancer types, ICB regimens, RECIST response labeling, and cBioPortal identifiers. |
| [`03_code_architecture_and_refactoring_guide.md`](./03_code_architecture_and_refactoring_guide.md) | Detailed code map tracing each dataset through downloading, preprocessing, integration, deconvolution, regression, and visualization, with actionable step-by-step refactoring guides for AI agents. |

---

## 3. Quick Dataset Reference Table

### A. Single-Cell Reference Datasets (Reference Signature Pool)

| Dataset Identifier | Cancer Type | Platform | Real Cell Depth | Clinical Labels? | Primary Code Files |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `SadeFeldman_Melanoma` (`GSE120575`) | Melanoma (pre/post ICB) | Smart-seq2 | 16,288 | Yes (Res/NR, Pre/Post) | `01_build_reference.py`, `04_analyze_milopy.py` |
| `JerbyArnon` (`GSE115978`) | Melanoma | Smart-seq2 | 7,186 | No | `01c_build_dataset_reference.py`, `01d_build_combined_references.py` |
| `Ma_Liver` (`GSE125449`) | Hepatocellular Carcinoma | 10x Chromium | 5,115 | No | `01c_build_dataset_reference.py`, `01d_build_combined_references.py` |
| `Yost_BCC` (`GSE123813`) | Basal Cell Carcinoma | 10x Chromium | 3,500 | No | `01c_build_dataset_reference.py`, `01d_build_combined_references.py` |
| `Maynard_NSCLC` (2020) | Non-Small Cell Lung | 10x Chromium | 3,000 | No | `01c_build_dataset_reference.py`, `01d_build_combined_references.py` |
| `Pelka_CRC` (`GSE178341`) | Colorectal Cancer | 10x Chromium | ~65,000 | No | `packages/singlecellrnasignature`, `_single_cell_datasets.py` |
| `Azizi_BRCA` | Breast Cancer | 10x Chromium | ~45,000 | No | `packages/singlecellrnasignature`, `_single_cell_datasets.py` |
| `Qian_PanCancer` | Pan-Cancer (Lung/CRC/OV/BRCA) | 10x Chromium | ~200,000 | No | `packages/singlecellrnasignature`, `_single_cell_datasets.py` |
| `Cheng_PanCancer` | Pan-Cancer T-Cells | 10x Chromium | ~390,000 | No | `packages/singlecellrnasignature`, `_single_cell_datasets.py` |
| `Leader_NSCLC` | Non-Small Cell Lung | 10x Chromium | ~35,000 | No | `packages/singlecellrnasignature`, `_single_cell_datasets.py` |
| `Kim_LUAD` | Lung Adenocarcinoma | 10x Chromium | ~40,000 | No | `packages/singlecellrnasignature`, `_single_cell_datasets.py` |
| `Becker_COAD` | Colorectal Cancer | 10x Chromium | ~30,000 | No | `packages/singlecellrnasignature`, `_single_cell_datasets.py` |
| `Khaliq_CC` | Colorectal Cancer | 10x Chromium | ~25,000 | No | `packages/singlecellrnasignature`, `_single_cell_datasets.py` |
| `Borcherding_ccRCC` | Clear Cell Renal Cell | 10x Chromium | ~25,000 | No | `packages/singlecellrnasignature`, `_single_cell_datasets.py` |
| `Sharma_HCC` | Liver Cancer | 10x Chromium | ~15,000 | No | `packages/singlecellrnasignature`, `_single_cell_datasets.py` |
| `Lu_HCC` | Hepatocellular Carcinoma | 10x Chromium | ~18,000 | No | `packages/singlecellrnasignature`, `_single_cell_datasets.py` |
| `Pu_PTC` | Papillary Thyroid | 10x Chromium | ~20,000 | No | `packages/singlecellrnasignature`, `_single_cell_datasets.py` |
| `Durante_UVM` | Uveal Melanoma | 10x Chromium | ~10,000 | No | `packages/singlecellrnasignature`, `_single_cell_datasets.py` |
| `Biermann_BrainMet` | Melanoma Brain Metastasis | 10x Chromium | ~12,000 | No | `packages/singlecellrnasignature`, `_single_cell_datasets.py` |
| `Vazquez_OV` | Ovarian Cancer | 10x Chromium | ~15,000 | No | `packages/singlecellrnasignature`, `_single_cell_datasets.py` |
| `Zhang_TNBC` (2021) | Triple-Negative Breast | 10x Chromium | ~18,000 | No | `packages/singlecellrnasignature`, `_single_cell_datasets.py` |
| `Zhang_Myeloid` (2022) | Pan-Cancer Myeloid | 10x Chromium | ~50,000 | No | `packages/singlecellrnasignature`, `_single_cell_datasets.py` |
| `GSE179994_Tcells` | Pan-Cancer T-Cells (ICI) | 10x Chromium | 150,849 | Partial (RDS meta) | `data/raw/GSE179994/` |

### B. Bulk Clinical Immunotherapy Validation Cohorts

| Cohort Identifier | Cancer Type | Treatment Regimen | Sample Size ($N$) | Responders / Non-Responders | Storage Path |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `Hugo-iAtlas` | Melanoma | anti-PD-1 (Pembrolizumab) | 27 | 14 / 13 | `scratch/lair/CBioPortalDataset-Hugo-iAtlas` |
| `Riaz-iAtlas` | Melanoma | anti-PD-1 (Nivolumab) | 107 | 34 / 73 | `scratch/lair/CBioPortalDataset-Riaz-iAtlas` |
| `Liu-iAtlas` | Melanoma | anti-PD-1 (Nivolumab / Pembro) | 122 | 48 / 74 | `scratch/lair/CBioPortalDataset-Liu-iAtlas` |
| `Gide-iAtlas` | Melanoma | anti-PD-1 +/- anti-CTLA-4 | 91 | 63 / 28 | `scratch/lair/CBioPortalDataset-Gide-iAtlas` |
| `Rosenberg-iAtlas` | Urothelial / Bladder | anti-PD-L1 (Atezolizumab) | 347 | 68 / 279 | `scratch/lair/CBioPortalDataset-Rosenberg-iAtlas` |
| `Padron-iAtlas` | Pancreatic Ductal | anti-PD-1 | 93 | 16 / 77 | `scratch/lair/CBioPortalDataset-Padron-iAtlas` |
| `Anders-iAtlas` | Breast (TNBC) | anti-PD-L1 | 31 | 11 / 20 | `scratch/lair/CBioPortalDataset-Anders-iAtlas` |
| `McDermott-iAtlas` | Renal Cell (ccRCC) | anti-PD-L1 +/- Bevacizumab | 263 | 97 / 166 | `scratch/lair/CBioPortalDataset-McDermott-iAtlas` |
| `Choueiri-iAtlas` | Renal Cell (ccRCC) | anti-PD-1 (Nivolumab) | 16 | 6 / 10 | `scratch/lair/CBioPortalDataset-Choueiri-iAtlas` |
| **Total** | **Pan-Cancer (5 Tissues)** | **All Regimens** | **1,097** | **357 / 740** | **9 Independent Studies** |
