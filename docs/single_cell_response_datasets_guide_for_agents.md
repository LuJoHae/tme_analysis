# Single-Cell Immunotherapy Response Datasets: Master Catalog & Query Guide for Agents

This document is the definitive reference manual for autonomous agents and downstream analysis pipelines querying single-cell RNA-sequencing (scRNA-seq / snRNA-seq) datasets with clinical immunotherapy response labels from the `tme_datasets` package.

---

## 1. Quick Start: Querying via `tme_datasets`

All datasets are registered in a centralized registry with deterministic paths, strict Ensembl Release 111 harmonization, and standardized observation columns.

### 1.1 Loading Any Response Cohort
```python
from returns.result import Success, Failure
from tme_datasets import load_dataset, get_dataset_metadata

# 1. Inspect metadata without loading heavy matrices
match get_dataset_metadata("GSE120575"):
    case Success(spec):
        print(f"Cohort: {spec.title}")
        print(f"Indication: {spec.cancer_type} | Platform: {spec.platform}")
        print(f"PMID: {spec.publication_pmid.value_or('N/A')}")
    case Failure(err):
        print(f"Metadata lookup failed: {err}")

# 2. Load standardized AnnData (sub-second cached H5AD read)
match load_dataset("GSE120575"):
    case Success(adata):
        print(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")
        print("Response breakdown:")
        print(adata.obs["clinical_response"].value_counts())
        print(f"Expression layers: {list(adata.layers.keys())}")
    case Failure(err):
        print(f"Load failed: {err}")
```

### 1.2 Declarative Discovery & Filtering
```python
from tme_datasets import list_datasets
from tme_datasets.models import Modality

# List all single-cell datasets with clinical response annotations
response_cohorts = list_datasets(modality=Modality.SINGLE_CELL, has_response=True)
for spec in response_cohorts:
    print(f"{spec.id:<22} | {spec.cancer_type:<15} | {spec.platform}")
```

### 1.3 Subsetting & In-Memory Filtering on Load
```python
# Filter on-the-fly without copying full matrices
res_responders = load_dataset(
    "GSE120575",
    subset={"clinical_response": "responder", "treatment_status": "pre-treatment"},
    subsample_n=5000,  # Deterministic reproducible downsampling
)
```

---

## 2. Standardized AnnData Contract & Observation Schema

Every AnnData object loaded via `tme_datasets` guarantees the following invariants:

```mermaid
classDiagram
    class StandardizedAnnData {
        +X: csr_matrix (Raw integer counts OR linear TPM if Smart-seq2)
        +layers['counts']: csr_matrix (Raw integer counts; omitted for TPM-only)
        +layers['tpm']: csr_matrix (Linear library-size scaled TPM, sum=10^6)
        +layers['log1p_norm']: csr_matrix (Natural log-transformed normalized: ln(1 + TPM))
        +obs['patient_id']: Categorical (Unique biological patient identifier)
        +obs['sample_id']: Categorical (Biospecimen / aliquot identifier)
        +obs['treatment_status']: Categorical (pre-treatment, on-treatment, post-treatment)
        +obs['clinical_response']: Categorical (responder, non-responder, stable, not-evaluable)
        +obs['clinical_response_raw']: Categorical (Verbatim ground-truth label from publication)
        +obs['cell_type']: Categorical (Curated or author cell lineage / cell state)
        +obs['total_counts']: float32 (Total UMI or count depth per cell)
        +obs['n_genes_by_counts']: int32 (Number of detected non-zero genes)
        +obs['pct_counts_mt']: float32 (Mitochondrial read percentage)
        +var_names: Index (HGNC Gene Symbol or Ensembl stable ID)
        +var['gene_name']: string (HGNC Gene Symbol)
        +var['gene_id']: string (Ensembl ID ENSG00000...)
        +uns['expression_type']: string ('counts' or 'tpm')
        +uns['is_smartseq2']: bool (True if full-length Smart-seq2)
        +uns['is_single_nucleus']: bool (True if snRNA-seq)
    }
```

### Response Columns Invariant
- **`adata.obs["clinical_response"]`**: Standardized binary/ternary classification:
  - `"responder"`: Patients achieving RECIST Complete Response (`CR`), Partial Response (`PR`), Pathological Complete Response (`pCR`), or Major Pathological Response (`MPR`).
  - `"non-responder"`: Patients with Progressive Disease (`PD`), Residual Disease (`RD`), Non-Major Pathological Response (`NMPR`), or clinical non-response.
  - `"stable"`: Patients with Stable Disease (`SD`) or intermediate response (e.g. Medium pathological response).
  - `"not-evaluable"`: Control samples, healthy donors, or unevaluated baseline timepoints.
- **`adata.obs["clinical_response_raw"]`**: Preserves the verbatim ground-truth string from the source table (`"High"`, `"MPR"`, `"pCR"`, `"CR: complete response"`, etc.).

---

## 3. Master Summary of the 9 Verified Clinical Response Cohorts

Following rigorous audit of 350 public single-cell cohorts across GEO and CELLxGENE, exactly **9 premier cohorts** possess verified patient-level or cell-level clinical immunotherapy response outcomes:

| Index | Accession / Canonical ID | Indication | Technology | Patients ($N$) | Cells | Standardized Response Breakdown | Primary Publication / DOI |
| :---: | :--- | :--- | :--- | :---: | :---: | :--- | :--- |
| **1** | **GSE120575** (`SadeFeldman`) | Melanoma | Smart-seq2 | 39 | 16,291 | 14 R / 25 NR | Sade-Feldman et al., *Cell* 2018; PMID: 30388456 |
| **2** | **CELLxGENE_7b20c613** (`7b20c613`) | Pan-Cancer | 10x 3' | 123 | 332,972 | 41 R / 82 NR | Gondal et al., *Sci Data* 2025; 10.1038/s41597-025-04381-6 |
| **3** | **CELLxGENE_05a8c945** (`05a8c945`) | Colorectal | 10x 3'/5' | 77 | 1,107,923 | 48 R / 29 NR | Marteau et al., *Cancer Cell* 2026; 10.1016/j.ccell.2025.12.003 |
| **4** | **CELLxGENE_6f9de485** (`6f9de485`) | Breast (TNBC) | 10x 3' | 82 | 424,114 | 45 R / 37 NR | Yan et al., *Nature* 2026; 10.1038/s41586-026-10469-9 |
| **5** | **GSE207422** (`Hu_NSCLC`) | NSCLC | BD Rhapsody | 14 | 91,997 | 4 R / 10 NR | Hu et al., *Genome Med* 2023; PMID: 36869384 |
| **6** | **GSE243013** (`Liu_NSCLC`) | NSCLC | 10x 3' | 242 | 1,254,749 | 130 R / 112 NR | Liu et al., *Cell* 2025; PMID: 40147443 |
| **7** | **GSE233203** (`Hong_NSCLC`) | NSCLC | 10x 3'/5' | 7 | 40,142 | 3 R / 4 NR | Hong et al., *MedComm* 2025; PMID: 40843133 |
| **8** | **GSE200996** (`Luoma_HNSCC`) | HNSCC | 10x 5' | 12 | 245,253 | 4 R / 8 NR | Luoma et al., *Cell* 2022; PMID: 35803260 |
| **9** | **GSE316195** (`Raufi_PDAC`) | PDAC | snRNA-seq | 17 | 83,651 | 15 R / 2 NR | Raufi et al., *Nat Commun* 2026; PMID: 42811017 |

---

## 4. In-Depth Cohort Dossiers

### 1. GSE120575 (Sade-Feldman et al. 2018 — Melanoma)
- **Aliases**: `Sade-Feldman`, `SadeFeldman`, `GSE120575`
- **Citation**: Sade-Feldman M, et al. *Defining T Cell States Associated with Response to Checkpoint Immunotherapy in Melanoma.* **Cell**. 2018 Nov 1;175(4):998-1013.e20. [PMID: 30388456](https://pubmed.ncbi.nlm.nih.gov/30388456/), [DOI: 10.1016/j.cell.2018.10.038](https://doi.org/10.1016/j.cell.2018.10.038).
- **Tumor Indication**: Metastatic Cutaneous Melanoma.
- **Therapy Regimen**: Anti-PD-1 monotherapy (nivolumab or pembrolizumab), anti-CTLA-4 monotherapy (ipilimumab), or concurrent ipilimumab + nivolumab.
- **Biopsy Sampling**: Longitudinal baseline (`Pre`, 32 biopsies) and on-treatment / post-progression (`Post`, 16 biopsies).
- **Sequencing & Isolation**: Plate-based full-length Smart-seq2 of FACS-isolated $\text{CD45}^+$ tumor-infiltrating immune cells.
- **Scale**: 16,291 cells across 37 patients (48 individual tumor biopsies).
- **Response Definition**:
  - `responder` ($N = 5{,}564$ cells): Patients achieving RECIST Complete Response, Partial Response, or durable Stable Disease ($\ge 6$ months).
  - `non-responder` ($N = 10{,}727$ cells): Patients with Progressive Disease.
- **Layer & Matrix Invariant**:
  - Full-length Smart-seq2 does not produce integer UMI counts; linear TPM is the genuine physical quantity.
  - Linear TPM is preserved in `.X` and `adata.layers["tpm"]`.
  - Natural log-transformed values $\ln(1 + \text{TPM})$ are stored in `adata.layers["log1p_norm"]`.
  - Artificial integer `layers["counts"]` is explicitly omitted.
  - Provenance tag `adata.uns["expression_type"] = "tpm"` and `adata.uns["is_smartseq2"] = True`.
- **Query Example**:
  ```python
  from tme_datasets import load_dataset
  adata = load_dataset("GSE120575").unwrap()
  tpm_matrix = adata.layers["tpm"]  # Linear TPM for deconvolution
  log_matrix = adata.layers["log1p_norm"]  # Log-transformed for differential expression
  ```

---

### 2. CELLxGENE_7b20c613 (Pan-Cancer ICB Meta-Atlas — Gondal et al. 2025)
- **Aliases**: `CELLxGENE_7b20c613`, `7b20c613`, `CELLxGENE_7b20c613_PanCancer`
- **Citation**: Gondal et al. *scICML: a pan-cancer single-cell multi-omics dataset for immune checkpoint inhibitor outcomes.* **Scientific Data**. 2025; [PMID: 39837887](https://pubmed.ncbi.nlm.nih.gov/39837887/), [DOI: 10.1038/s41597-025-04381-6](https://doi.org/10.1038/s41597-025-04381-6).
- **Tumor Indication**: Pan-Cancer ICB Meta-Atlas (58.5% Breast, 23.0% Basal Cell Carcinoma, 8.7% ccRCC, 7.1% Melanoma, 2.6% Liver).
- **Therapy Regimen**: Anti-PD-1 (pembrolizumab, nivolumab) monotherapy and combination ICB.
- **Sequencing & Isolation**: 10x Genomics Chromium 3' (v2/v3) droplet scRNA-seq.
- **Scale**: 332,972 cells across 123 clinically annotated immunotherapy-treated patients (41 Responders, 82 Non-Responders).
- **Response Definition**:
  - `responder` ($N = 41$ patients): Harmonized favorable response (`CR`, `PR`, `R`, `Favourable`).
  - `non-responder` ($N = 82$ patients): Harmonized unfavorable response (`PD`, `NR`, `Unfavourable`).
  - Verbatim keys: `adata.obs["clinical_response_raw"]` stores original `Combined_outcome` and `outcome`.
- **Patient Identifier Safeguard**: Patient identifiers are prefixed with study accession in `adata.obs["patient_id"]` to guarantee zero patient-key collisions across integrated cohorts.
- **Query Example**:
  ```python
  adata = load_dataset("CELLxGENE_7b20c613").unwrap()
  t_cells = adata[adata.obs["cell_type"].str.contains("T cell")].copy()
  ```

---

### 3. CELLxGENE_05a8c945 (Colorectal Cancer Global Core Atlas — Marteau et al. 2026)
- **Aliases**: `CELLxGENE_05a8c945`, `05a8c945`, `CELLxGENE_05a8c945_CRC`
- **Citation**: Marteau P, et al. *Single-cell integration and multi-modal profiling reveals phenotypes of critical colorectal cancer continuum.* **Cancer Cell**. 2026; [PMID: 41453298](https://pubmed.ncbi.nlm.nih.gov/41453298/), [DOI: 10.1016/j.ccell.2025.12.003](https://doi.org/10.1016/j.ccell.2025.12.003). (CELLxGENE dataset accession `05a8c945-bc12-414f-960d-a31943bbcdd1`).
- **Tumor Indication**: Colorectal Carcinoma (both dMMR/MSI-H and pMMR/MSS).
- **Therapy Regimen**: Neoadjuvant pembrolizumab / nivolumab (anti-PD-1 monotherapy or anti-PD-1 + anti-CTLA-4 combination).
- **Sequencing & Isolation**: 10x Chromium 3' and 5' multi-omics.
- **Scale**: 1,107,923 cells across 77 patients with clinical response annotations (48 Responders, 29 Non-Responders).
- **Response Definition**:
  - `responder` ($N = 48$ patients): RECIST Complete Response (`CR: complete response`) or Partial Response (`PR: partial response`).
  - `non-responder` ($N = 29$ patients): Progressive Disease (`PD: progressive disease`).
  - Stable / non-evaluable: `SD: stable disease`, `NE: inevaluable`.
- **Memory Safety**: Filtered at load time via backed streaming to guarantee sub-5 GB RAM footprint despite the source atlas exceeding 30 GB on disk.
- **Query Example**:
  ```python
  adata = load_dataset("CELLxGENE_05a8c945").unwrap()
  # obs columns include: 'RECIST', 'treatment_response', 'treatment_status_before_resection'
  ```

---

### 4. CELLxGENE_6f9de485 (Triple-Negative Breast Cancer ARTEMIS Trial — Yan et al. 2026)
- **Aliases**: `CELLxGENE_6f9de485`, `6f9de485`, `CELLxGENE_6f9de485_Breast`
- **Citation**: Yan Y, et al. (Navin laboratory). *Single-cell and spatial dissection of neoadjuvant chemotherapy response in triple-negative breast cancer (ARTEMIS trial).* **Nature**. 2026; [PMID: 41586026](https://pubmed.ncbi.nlm.nih.gov/41586026/), [DOI: 10.1038/s41586-026-10469-9](https://doi.org/10.1038/s41586-026-10469-9).
- **Tumor Indication**: Triple-Negative Breast Cancer (TNBC).
- **Therapy Regimen**: Neoadjuvant Chemotherapy (AC: Doxorubicin + Cyclophosphamide, ARTEMIS trial).
- **Sequencing & Isolation**: 10x Genomics Chromium 3' v3.
- **Scale**: 424,114 cells across 82 patients (45 Responders / pCR, 37 Non-Responders / RD).
- **Response Definition**:
  - `responder` ($N = 45$ patients): Pathological Complete Response (`pCR`).
  - `non-responder` ($N = 37$ patients): Residual Disease (`RD`).
- **Query Example**:
  ```python
  adata = load_dataset("CELLxGENE_6f9de485").unwrap()
  myeloid = adata[adata.obs["cell_type"].isin(["Macrophage", "Dendritic cell"])].copy()
  ```

---

### 5. GSE207422 (NSCLC Neoadjuvant Tislelizumab/Camrelizumab Trial — Hu et al. 2023)
- **Aliases**: `GSE207422`, `Hu_NSCLC`, `GSE207422_NSCLC`
- **Citation**: Hu J, et al. *Dissecting theimpact of neoadjuvant PD-1 blockade and chemotherapy on the tumor microenvironment of non-small cell lung cancer.* **Genome Medicine**. 2023; [PMID: 36869384](https://pubmed.ncbi.nlm.nih.gov/36869384/), [DOI: 10.1186/s13073-023-01164-9](https://doi.org/10.1186/s13073-023-01164-9).
- **Tumor Indication**: Resectable Non-Small Cell Lung Cancer (NSCLC).
- **Therapy Regimen**: Neoadjuvant anti-PD-1 (tislelizumab / camrelizumab) + platinum-doublet chemotherapy.
- **Sequencing & Isolation**: BD Rhapsody Targeted & Whole-Transcriptome single-cell RNA-seq.
- **Scale**: 91,997 cells across 14 patients (4 Responders, 10 Non-Responders).
- **Response Definition**:
  - `responder` ($N = 4$ patients): Major Pathological Response (`MPR`, $\le 10\%$ viable residual tumor cells) and Complete Pathological Response (`pCR`).
  - `non-responder` ($N = 10$ patients): Non-Major Pathological Response (`NMPR`, $>10\%$ viable tumor).
  - Radiographic RECIST: `PR`, `SD` documented in `adata.obs["clinical_response_raw"]`.
- **Query Example**:
  ```python
  adata = load_dataset("GSE207422").unwrap()
  ```

---

### 6. GSE243013 (NSCLC Megascale Multicenter Trial — Liu et al. 2025)
- **Aliases**: `GSE243013`, `Liu_NSCLC`, `GSE243013_NSCLC`
- **Citation**: Liu Z, et al. *Single-cell atlas of non-small cell lung cancer treated with neoadjuvant immunochemotherapy.* **Cell**. 2025; [PMID: 40147443](https://pubmed.ncbi.nlm.nih.gov/40147443/), [DOI: 10.1016/j.cell.2025.03.018](https://doi.org/10.1016/j.cell.2025.03.018).
- **Tumor Indication**: Non-Small Cell Lung Cancer (NSCLC).
- **Therapy Regimen**: Neoadjuvant anti-PD-(L)1 (sintilimab / tislelizumab / atezolizumab) + chemotherapy.
- **Sequencing & Isolation**: 10x Chromium 3' v3 multicenter high-throughput scRNA-seq.
- **Scale**: **1,254,749 cells across 242 patients** (130 Responders, 112 Non-Responders; largest single-cohort clinical ICB single-cell resource published).
- **Response Definition**:
  - `responder` ($N = 130$ patients): Major Pathological Response (`MPR`).
  - `non-responder` ($N = 112$ patients): Non-Major Pathological Response (`non-MPR`).
- **Query Example**:
  ```python
  adata = load_dataset("GSE243013").unwrap()
  print(f"Total NSCLC cells: {adata.n_obs} across {adata.obs['patient_id'].nunique()} patients")
  ```

---

### 7. GSE233203 (NSCLC ABCP Quad Therapy Trial — Hong et al. 2025)
- **Aliases**: `GSE233203`, `Hong_NSCLC`, `GSE233203_NSCLC`
- **Citation**: Hong W, et al. *Single-cell dynamics under atezolizumab, bevacizumab, carboplatin, and paclitaxel in advanced non-small cell lung cancer.* **MedComm**. 2025; [PMID: 40843133](https://pubmed.ncbi.nlm.nih.gov/40843133/), [DOI: 10.1002/mco2.70330](https://doi.org/10.1002/mco2.70330).
- **Tumor Indication**: Advanced / Metastatic NSCLC.
- **Therapy Regimen**: Atezolizumab + Bevacizumab + Carboplatin + Paclitaxel (ABCP quad therapy).
- **Sequencing & Isolation**: 10x Chromium 5' and 3' single-cell RNA-seq.
- **Scale**: 40,142 single cells across 7 patients (3 Responders, 4 Non-Responders).
- **Response Definition**:
  - `responder` ($N = 3$ patients): Confirmed RECIST Objective Response (`PR`).
  - `non-responder` ($N = 4$ patients): Stable Disease (`SD`) or Progressive Disease (`PD`).
- **Query Example**:
  ```python
  adata = load_dataset("GSE233203").unwrap()
  ```

---

### 8. GSE200996 (HNSCC TIL Atlas — Luoma et al. 2022)
- **Aliases**: `GSE200996`, `Luoma_HNSCC`, `GSE200996_HNSCC`
- **Citation**: Luoma AM, et al. *Tissue-resident memory and circulating T cells are early responders to pre-surgical cancer immunotherapy.* **Cell**. 2022; [PMID: 35803260](https://pubmed.ncbi.nlm.nih.gov/35803260/), [DOI: 10.1016/j.cell.2022.06.018](https://doi.org/10.1016/j.cell.2022.06.018).
- **Tumor Indication**: Oral Cavity & Laryngeal Head and Neck Squamous Cell Carcinoma (HNSCC).
- **Therapy Regimen**: Neoadjuvant anti-PD-1 (nivolumab) monotherapy vs. nivolumab + ipilimumab combination.
- **Sequencing & Isolation**: 10x Genomics 5' immune profiling with paired TCR $\alpha\beta$ sequencing.
- **Scale**: 245,253 cells across 12 neoadjuvant ICB patients with response evaluations (4 Responders, 8 Non-Responders).
- **Response Definition**:
  - `responder` ($N = 4$ patients): High pathological treatment response (`Path_response: High`, $>50\%$ tumor necrosis/regression).
  - `non-responder` ($N = 8$ patients): Medium or Low pathological response (`Path_response: Medium/Low`, $<50\%$).
- **Optimization Note**: Raw GEO archive contains 30.2M unfiltered droplets across 204 `.h5` files. `tme_datasets` automatically performs droplet thresholding ($\ge 200$ counts) and case-insensitive barcode matching against author metadata tables.
- **Query Example**:
  ```python
  adata = load_dataset("GSE200996").unwrap()
  cd8_trm = adata[adata.obs["clinical_response"].isin(["responder", "non-responder"])].copy()
  ```

---

### 9. GSE316195 (PDAC Combination Trial — Raufi et al. 2026)
- **Aliases**: `GSE316195`, `Raufi_PDAC`, `GSE316195_PDAC`
- **Citation**: Raufi AG, et al. (Izar and Manji laboratories). *Single-nucleus transcriptomics of human pancreatic ductal adenocarcinoma undergoing CXCR4 antagonist and anti-PD-1 therapy.* **Nature Communications**. 2026; [PMID: 42811017](https://pubmed.ncbi.nlm.nih.gov/42811017/), [DOI: 10.1038/s41467-026-69745-z](https://doi.org/10.1038/s41467-026-69745-z).
- **Tumor Indication**: Metastatic Pancreatic Ductal Adenocarcinoma (PDAC).
- **Therapy Regimen**: Motixafortide (CXCR4 antagonist) + cemiplimab (anti-PD-1) + gemcitabine / nab-paclitaxel chemotherapy.
- **Sequencing & Isolation**: 10x Chromium Single-Nucleus RNA-seq (`snRNA-seq`) from flash-frozen core needle liver metastases.
- **Scale**: 83,651 single nuclei across 17 patients (15 Responders, 2 Non-Responders).
- **Response Definition**:
  - `responder` ($N = 15$ patients): Confirmed RECIST Partial Response (`PR`).
  - `non-responder` ($N = 2$ patients): Progressive Disease (`PD`) or Stable Disease (`SD`).
  - Verbatim response: `PR`, `SD`, `PD` in `adata.obs["clinical_response_raw"]`.
- **QC Parameters**: Tagged with `adata.uns["is_single_nucleus"] = True` (mitochondrial cutoff automatically tightened to $10.0\%$).
- **Query Example**:
  ```python
  adata = load_dataset("GSE316195").unwrap()
  assert adata.uns["is_single_nucleus"] is True
  ```

---

## 5. Overview of Tier 1 ICB-Treated Cohorts (Class C)

The following cohorts represent bona fide human clinical trials where patients were treated with immune checkpoint blockade, but public records provide **overall cohort treatment setting, longitudinal timepoints, or survival, rather than discrete per-cell binary response labels**:

| Accession | Indication | Therapy | Why Not Tier 0 (Class C Rationale) |
| :--- | :--- | :--- | :--- |
| **GSE246613** | Breast (TNBC) | Pembrolizumab + Radiotherapy | Metadata provides treatment cycle (`cycle_1`, `cycle_2`) and trajectory clusters, but individual binary RECIST labels are withheld in GEO. |
| **GSE302453** | Breast | Neoadjuvant anti-PD-1 | Cohort is ICB-treated, but published metadata matrix lacks individual patient response column. |
| **GSE270464** | Melanoma | Anti-PD-1 | Acral/mucosal melanoma trial; individual patient outcome classifications not provided in supplementary tables. |
| **GSE276139** | NSCLC | Anti-PD-1 + Chemotherapy | CSF cells from leptomeningeal metastases; treated cohort without RECIST stratification. |
| **GSE220313** | ccRCC | Anti-PD-1 (Ex Vivo / In Vivo) | In vitro stimulation of sorted $\text{CD45}^+$ cells from patients on nivolumab; no clinical response table. |
| **GSE254498** | ccRCC | Anti-PD-1 | VHL disease cohort; descriptive case series without RECIST responder stratification. |

To query Tier 1 cohorts for microenvironmental or atlas studies:
```python
tier1_cohorts = list_datasets(modality=Modality.SINGLE_CELL, tier="Tier 1")
```

---

## 6. Guidelines for Downstream Agent Workflows

### 1. Training Response Predictors / Biomarkers
- **Strictly partition by patient**: Never split cells randomly across train and test splits. Cross-validation MUST group by `adata.obs["patient_id"]` to avoid severe data leakage:
  ```python
  from sklearn.model_selection import GroupKFold
  gkf = GroupKFold(n_splits=5)
  groups = adata.obs["patient_id"].values
  ```
- **Filter to Pre-Treatment**: When predicting baseline therapeutic response, filter out on-treatment biopsies:
  ```python
  baseline_adata = adata[adata.obs["treatment_status"].isin(["baseline", "pre-treatment"])].copy()
  ```

### 2. Cell-Type Specific Pseudobulking
- Aggregate single cells into patient-level pseudobulk profiles using `tme_datasets.harmonization`:
  ```python
  from tme_datasets.harmonization.pseudobulk import compute_cell_type_pseudobulk
  pbulk_df = compute_cell_type_pseudobulk(adata, cell_type_col="cell_type", patient_col="patient_id")
  ```

### 3. Deconvolution Signature Construction
- For signature reference matrix $\Phi$, use linear TPM from `adata.layers["tpm"]` or `adata.X`:
  ```python
  # Ensure linear space (not log-transformed)
  lin_expr = adata.layers["tpm"] if "tpm" in adata.layers else adata.X
  ```
