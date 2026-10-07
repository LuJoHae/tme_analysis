# Single-Cell Reference Datasets Catalog

> [!IMPORTANT]
> For agents seeking the complete list and programmatic guide for **single-cell datasets with clinical immunotherapy response labels** queried via `tme_datasets`, refer to:
> [`docs/single_cell_response_datasets_guide_for_agents.md`](file:///Users/halu/Code/tme_analysis/docs/single_cell_response_datasets_guide_for_agents.md).

This document provides extensive details on all single-cell RNA-sequencing reference datasets in the repository, both those with clinical response labels and those without.

---

## 1. Single-Cell Response Labels vs. Reference Construction

In single-cell deconvolution benchmarking:
- **Reference Signature Generation ($\Phi$)**: Does **NOT** require response labels. A single-cell dataset is used purely to define immune cell states (via graph clustering) and compute linear expression profiles:
  $$\Phi_{k,g} = \frac{\bar{x}_{k,g}}{\sum_{g'} \bar{x}_{k,g'}}$$
  Any high-quality single-cell dataset covering the tumor microenvironment (T-cells, B-cells, NK, macrophages, DCs, endothelial, fibroblasts) can serve as a reference.
- **Bulk Clinical Testing**: Clinical response labels (Responder vs. Non-Responder) are evaluated on the bulk RNA-seq cohorts ($N=1,097$), not the single-cell reference.
- **Single-cell Milo Differential Abundance (Step 4)**: The only single-cell dataset where response labels are actively utilized is **Sade-Feldman et al. (`GSE120575`)**, where cells have annotated patient response (`Responder` vs `Non-Responder`) and biopsy timepoint (`Pre` vs `Post`). This provides an orthogonal, single-cell-only ground truth for comparing bulk-inferred odds ratios ($\beta_{\text{bulk}}$) against single-cell log2 fold-changes ($\text{LFC}_{\text{Milo}}$).

---

## 2. Locally Unpacked Primary Single-Cell Datasets

These 5 datasets are downloaded, unpacked, and fully ready in the repository workspace without any artificial subsampling:

### A. Sade-Feldman et al. 2018 (`GSE120575`)
* **Study**: *Defining T Cell States Associated with Response to Checkpoint Immunotherapy in Melanoma* (*Cell*, 2018).
* **Tissue**: Metastatic Melanoma (32 patients, 48 tumor biopsies).
* **Sequencing Platform**: Full-length Smart-seq2 (plate-based).
* **Cell Count**: **16,288 single cells** (all immune cells: CD45+ sorted).
* **Clinical Labels**:
  - Response: Responder (CR/PR, long-term SD) vs. Non-Responder (PD).
  - Treatment: Anti-PD-1 monotherapy, anti-CTLA-4 monotherapy, or combination.
  - Timepoint: Baseline (`Pre`) and on-treatment / post-progression (`Post`).
* **Disk Locations**:
  - Raw files: `data/raw/GSE120575/GSE120575_tpm.txt.gz`, `data/raw/GSE120575/GSE120575_meta.txt.gz`
  - Parquet format: `scratch/GSE120575/gse120575_tpm.parquet`, `scratch/GSE120575/gse120575_tpm_cell_metadata.parquet`
  - Preprocessed reference: `output/sade_feldman_deconv_validation/reference_phi.parquet`
* **Code Implementation**:
  - Load script: [`scripts/sade_feldman_deconv_validation/01_build_reference.py`](file:///Users/halu/Code/tme_analysis/scripts/sade_feldman_deconv_validation/01_build_reference.py#L40-L100)
  - Processing: Confounding genes filtered out (ribosomal, mitochondrial, HLA/immunoglobulin variable chains, melanin/tumor markers), linear TPM converted to cluster means, normalized to unit simplex.

### B. Jerby-Arnon et al. 2018 (`GSE115978`)
* **Study**: *A Cancer Cell Program Promotes T Cell Exclusion and Resistance to Checkpoint Blockade* (*Cell*, 2018).
* **Tissue**: Metastatic Melanoma (33 patients).
* **Sequencing Platform**: Smart-seq2.
* **Cell Count**: **7,186 single cells** (immune, stromal, and malignant).
* **Response Labels**: None (pre-treatment melanoma cross-sectional cohort).
* **Disk Location**:
  - Raw files: `data/raw/GSE115978/GSE115978_meta.gz`, `data/raw/GSE115978/GSE115978_tpm.gz`
* **Code Implementation**:
  - Ingestion: [`scripts/sade_feldman_deconv_validation/01c_build_dataset_reference.py`](file:///Users/halu/Code/tme_analysis/scripts/sade_feldman_deconv_validation/01c_build_dataset_reference.py#L140-L192)
  - Function: `load_jerby_arnon(raw_dir: Path)`
  - Processing: Reads CSV linear TPM, filters confounding gene families, constructs AnnData, runs log1p, PCA (30 PCs), k-NN (k=15), and Leiden clustering across 8 resolutions (`0.25` to `2.00`).

### C. Ma et al. 2019 (`GSE125449`)
* **Study**: *Single-Cell RNA-seq Reveals Heterogeneity of Tumor Microenvironment in Liver Cancer* (*Cell*, 2019).
* **Tissue**: Hepatocellular Carcinoma (HCC) and Intrahepatic Cholangiocarcinoma (ICC).
* **Sequencing Platform**: 10x Genomics Chromium (v2).
* **Cell Count**: **5,115 single cells** (tumor-infiltrating immune and stromal cells).
* **Response Labels**: None (surgical resection, treatment-naive / untreated).
* **Disk Location**:
  - Raw files: `data/raw/GSE125449/GSE125449_set1_matrix.gz`, `GSE125449_set1_genes.gz`, `GSE125449_set1_barcodes.gz`
* **Code Implementation**:
  - Ingestion: [`scripts/sade_feldman_deconv_validation/01c_build_dataset_reference.py`](file:///Users/halu/Code/tme_analysis/scripts/sade_feldman_deconv_validation/01c_build_dataset_reference.py#L226-L280)
  - Function: `load_ma_liver(raw_dir: Path)`
  - Processing: Reads MTX sparse matrix, normalizes to linear CPM ($10^6$ library depth scaling), constructs AnnData, clusters across resolutions 0.25–2.00.

### D. Yost et al. 2019 (`GSE123813` / `GSE123139`)
* **Study**: *Clonal replacement of tumor-specific T cells following PD-1 blockade* (*Nature Medicine*, 2019).
* **Tissue**: Basal Cell Carcinoma (BCC) and Squamous Cell Carcinoma (SCC).
* **Sequencing Platform**: 10x Genomics Chromium (5' immune profiling).
* **Cell Count**: **3,500 single cells** (BCC CD45+ immune infiltrate subset; >15,000 total across BCC + SCC).
* **Response Labels**: None in standard release (contains pre/post biopsy pairs).
* **Disk Location**:
  - Raw files: `data/raw/GSE123813/GSE123813_bcc_counts.txt.gz`, `GSE123813_scc_counts.txt.gz`
  - Series matrix: `data/raw/GSE123139/GSE123139_matrix.txt.gz`
* **Code Implementation**:
  - Ingestion: [`scripts/sade_feldman_deconv_validation/01c_build_dataset_reference.py`](file:///Users/halu/Code/tme_analysis/scripts/sade_feldman_deconv_validation/01c_build_dataset_reference.py#L282-L340)
  - Function: `load_yost(raw_dir: Path)`
  - Processing: Reads UMI count matrix, scales to linear CPM, filters confounding genes, clusters across 8 resolutions.

### E. Maynard et al. 2020 (NSCLC)
* **Study**: *Therapy-Induced Evolution of Human Lung Cancer Revealed by Single-Cell RNA-Seq* (*Cell*, 2020).
* **Tissue**: Advanced Non-Small Cell Lung Cancer (targeted therapy / ICB longitudinal cohort across 45 patients).
* **Sequencing Platform**: Smart-seq2 (plate-based full-length single-cell RNA-seq).
* **Cell Count**: **27,489 single cells** across 49 biopsies (67.1M non-zero measurements).
* **Response Labels**: Longitudinal biopsy timepoints (`biopsy_time_status`: treatment-naive TN, residual disease RD, progressive disease PD), annotated response status, and patient metadata.
* **Disk Location**:
  - H5AD file: `data/manual_download/Maynard_NSCLC.h5ad` (or preprocessed cache `data/preprocessed/Maynard_NSCLC.h5ad`)
  - Raw source files: `data/raw/Maynard_NSCLC/S01_datafinal.csv` (1.47 GB), `data/raw/Maynard_NSCLC/S01_metacells.csv` (10.6 MB)
* **Code Implementation**:
  - Ingestion & Builder: [`packages/tme_datasets/src/tme_datasets/providers/single_cell.py`](file:///Users/halu/Code/tme_analysis/packages/tme_datasets/src/tme_datasets/providers/single_cell.py)
  - Function: `download_and_build_maynard_full(raw_dir, output_h5ad)` / `load_maynard(repo_root)`
  - Processing: Memory-efficient streaming parsing into sparse CSR AnnData, standardized clinical metadata (`patient`, `sample`, `timepoint`, `response_binary`), and TME lineage scoring.

### F. Tietscher et al. (`GSE179994`)
* **Study**: Pan-cancer tumor-infiltrating T-cell atlas across multiple cancer types.
* **Sequencing Platform**: 10x Genomics Chromium.
* **Cell Count**: **150,849 single T cells**.
* **Disk Location**:
  - Raw files: `data/raw/GSE179994/GSE179994_counts.gz` (gzip-compressed R RDS matrix, 441 MB gzip, 2.4 GB uncompressed), `GSE179994_meta.gz` (metadata with cellid, patient, sample, celltype, cluster).

---

## 3. Pan-Cancer Atlas Cohorts (`singlecellrnasignature`)

The repository contains an integrated pipeline that draws from 17 single-cell datasets defined in `packages/singlecellrnasignature` and `packages/single_cell_datasets`. None of these datasets have response labels; they are cross-sectional tumor microenvironment atlases.

| Dataset Identifier | Study Reference | Tumor Type | Tech | Full Depth (Est.) | Download Class |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `PelkaSpatiallyOrganizedMulticellular2021` | Pelka et al. *Cell* 2021 | Colorectal Cancer | 10x v3 | ~65,000 | `singlecellrnasignature.adata.PelkaSpatiallyOrganizedMulticellular2021Adata` |
| `AziziSingleCellMapDiverse2018` | Azizi et al. *Cell* 2018 | Breast Cancer | 10x v2 | ~45,000 | `singlecellrnasignature.adata.AziziSingleCellMapDiverse2018Adata` |
| `QianPancancerBlueprintHeterogeneous2020` | Qian et al. *Cell Res* 2020 | Pan-Cancer (Lung, CRC, Ovary, Breast) | 10x v2/v3 | ~200,000 | `singlecellrnasignature.adata.QianPancancerBlueprintHeterogeneous2020aAdata` |
| `ChengPancancerSinglecellTranscriptional2021` | Cheng et al. *Cell* 2021 | Pan-Cancer T-Cells | 10x v3 | ~390,000 | `singlecellrnasignature.adata.ChengPancancerSinglecellTranscriptional2021Adata` |
| `LeaderSinglecellAnalysisHuman2021` | Leader et al. *Cancer Cell* 2021 | Non-Small Cell Lung | 10x v3 | ~35,000 | `singlecellrnasignature.adata.LeaderSinglecellAnalysisHuman2021Adata` |
| `KimSinglecellRNASequencing2020` | Kim et al. *Nat Commun* 2020 | Lung Adenocarcinoma | 10x v2 | ~40,000 | `singlecellrnasignature.adata.KimSinglecellRNASequencing2020Adata` |
| `BeckerSinglecellAnalysesDefine2022` | Becker et al. *Cell* 2022 | Colorectal Cancer | 10x v3 | ~30,000 | `singlecellrnasignature.adata.BeckerSinglecellAnalysesDefine2022Adata` |
| `KhaliqRefiningColorectalCancer2022` | Khaliq et al. *iScience* 2022 | Colorectal Cancer | 10x v3 | ~25,000 | `singlecellrnasignature.adata.KhaliqRefiningColorectalCancer2022Adata` |
| `BorcherdingMappingImmuneEnvironment2021` | Borcherding et al. *Cancer Discov* 2021 | Clear Cell Renal Cell | 10x v3 | ~25,000 | `singlecellrnasignature.adata.BorcherdingMappingImmuneEnvironment2021Adata` |
| `SharmaOncofetalReprogrammingEndothelial2020` | Sharma et al. *Cell* 2020 | Hepatocellular Carcinoma | 10x v2 | ~15,000 | `singlecellrnasignature.adata.SharmaOncofetalReprogrammingEndothelial2020Adata` |
| `LuSinglecellAtlasMulticellular2022` | Lu et al. *Gut* 2022 | Liver Cancer | 10x v3 | ~18,000 | `singlecellrnasignature.adata.LuSinglecellAtlasMulticellular2022Adata` |
| `PuSinglecellTranscriptomicAnalysis2021` | Pu et al. *Genomics Proteomics Bioinform* 2021 | Papillary Thyroid | 10x v3 | ~20,000 | `singlecellrnasignature.adata.PuSinglecellTranscriptomicAnalysis2021Adata` |
| `DuranteSinglecellAnalysisReveals2020` | Durante et al. *Nat Commun* 2020 | Uveal Melanoma | 10x v2 | ~10,000 | `singlecellrnasignature.adata.DuranteSinglecellAnalysisReveals2020Adata` |
| `BiermannDissectingTreatmentnaiveEcosystem2022` | Biermann et al. *Cancer Cell* 2022 | Melanoma Brain Met | 10x v3 | ~12,000 | `singlecellrnasignature.adata.BiermannDissectingTreatmentnaiveEcosystem2022Adata` |
| `VazquezOvarianCancerMutational2022` | Vazquez et al. *Cell Rep* 2022 | Ovarian Cancer | 10x v3 | ~15,000 | `singlecellrnasignature.adata.VazquezOvarianCancerMutational2022Adata` |
| `ZhangSinglecellAnalysesReveal2021` | Zhang et al. *Cancer Cell* 2021 | Triple-Negative Breast | 10x v3 | ~18,000 | `singlecellrnasignature.adata.ZhangSinglecellAnalysesReveal2021Adata` |
| `ZhangSinglecellAnalysisReveals2022` | Zhang et al. *Cell* 2022 | Pan-Cancer Myeloid | 10x v3 | ~50,000 | `singlecellrnasignature.adata.ZhangSinglecellAnalysisReveals2022Adata` |

---

## 4. Origin of the "2,083 Cells" Subsampling Cap

In the pre-compiled `Combined-Atlas` (41,284 cells), cohorts 2 through 10 in Table 1 each had exactly 2,083 cells.

### Mathematical Origin:
In [`scripts/sade_feldman_deconv_validation/01b_build_integrated_reference.py`](file:///Users/halu/Code/tme_analysis/scripts/sade_feldman_deconv_validation/01b_build_integrated_reference.py#L37-L85):
```python
n_subsample_atlas: int = 25000

# Stratified sampling across 12 cancer codes
c_codes = obs["cancer_code"].unique()  # exactly 12 cancer types
n_per_code = max(500, n_sample // len(c_codes))  # 25,000 // 12 = 2,083

for code in c_codes:
    sub_indices = np.where(obs["cancer_code"] == code)[0]
    if len(sub_indices) > n_per_code:
        chosen = np.random.choice(sub_indices, size=n_per_code, replace=False)  # draws 2,083 cells!
```
- $12 \times 2{,}083 = 24{,}996$ cells sampled from the 10x atlas.
- $24{,}996 + 16{,}288 \text{ (Sade-Feldman)} = \mathbf{41{,}284}\text{ cells}$.

### How to Bypass This Cap for Full Datasets:
When another agent refactors reference integration:
1. Remove or set `n_sample = None` in `01b_build_integrated_reference.py`.
2. In `01d_build_combined_references.py`, add `sf` (Sade-Feldman) to the component dataset loader (`load_sade_feldman`) so that all 5 local datasets are integrated without any cap:
   $$16{,}288 + 7{,}186 + 5{,}115 + 3{,}500 + 3{,}000 = \mathbf{35{,}089}\text{ real single cells}.$$
