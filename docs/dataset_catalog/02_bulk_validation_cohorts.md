# Bulk Clinical Validation Cohorts Catalog

This document details the bulk clinical transcriptomic cohorts used to evaluate whether inferred cell-state abundances predict patient immunotherapy response.

---

## 1. Role in the Deconvolution & Predictive Pipeline

* **Input to Deconvolution**: Bulk tumor RNA-sequencing expression matrices $X \in \mathbb{R}^{G \times S}$ (genes $\times$ tumor samples).
* **Deconvolution Method**: BayesPrism / InstaPrism non-negative least squares (NNLS), decomposing bulk mixtures into cell state fractions $\theta_{s,k} \in [0, 1]$ where $\sum_{k=1}^K \theta_{s,k} = 1$.
* **Response Labels**: Every sample in these cohorts has annotated clinical outcome data (RECIST response: CR, PR, SD, PD; binarized to `Responder` vs. `Non-Responder`).
* **Statistical Modeling (Step 3)**:
  - Univariate logistic regression: $\text{logit}(P(\text{Response}=1)) = \alpha + \beta_k \cdot \theta_k$.
  - Multivariate elastic-net / L2 penalized logistic regression combining all $K$ cell fractions.
  - Stratified ROC-AUC calculation across **Melanoma**, **Pan-Cancer**, and individual cohorts.

---

## 2. Primary 9 iAtlas / cBioPortal Validation Cohorts ($N = 1,097$)

All 9 cohorts are stored and cached locally in `scratch/lair/CBioPortalDataset-<Cohort>-iAtlas`:

### A. Hugo et al. 2016 (`Hugo-iAtlas`)
* **Publication**: *Genomic and Transcriptomic Features of Response to Anti-PD-1 Therapy in Metastatic Melanoma* (*Cell*, 2016; DOI: [10.1016/j.cell.2016.02.065](https://doi.org/10.1016/j.cell.2016.02.065)).
* **GEO Accession**: `GSE78220` | **cBioPortal Identifier**: `skcm_hugo_2016`
* **Cancer Type**: Metastatic Cutaneous Melanoma.
* **Treatment**: Anti-PD-1 monotherapy (Pembrolizumab).
* **Sample Size**: **$N = 27$** pre-treatment tumor biopsies.
* **Response Breakdown**: 14 Responders (51.9%) / 13 Non-Responders (48.1%).
* **Disk Location**: `scratch/lair/CBioPortalDataset-Hugo-iAtlas`
  - Expression: `data_mrna_seq_fpkm.txt` (or RSEM TPM)
  - Clinical: `data_clinical_sample.txt`, `data_clinical_patient.txt`

### B. Riaz et al. 2017 (`Riaz-iAtlas`)
* **Publication**: *Tumor and Microenvironment Evolution during Immunotherapy with Nivolumab* (*Cell*, 2017; DOI: [10.1016/j.cell.2017.09.028](https://doi.org/10.1016/j.cell.2017.09.028)).
* **GEO Accession**: `GSE91061` | **cBioPortal Identifier**: `skcm_riaz_2017`
* **Cancer Type**: Advanced Melanoma (Ipilimumab-naive and Ipilimumab-pretreated).
* **Treatment**: Anti-PD-1 monotherapy (Nivolumab).
* **Sample Size**: **$N = 107$** samples (includes pre-treatment and early on-treatment).
* **Response Breakdown**: 34 Responders (31.8%) / 73 Non-Responders (68.2%).
* **Disk Location**: `scratch/lair/CBioPortalDataset-Riaz-iAtlas`

### C. Liu et al. 2019 (`Liu-iAtlas`)
* **Publication**: *Integrative molecular and clinical modeling of clinical outcomes to PD1 blockade in patients with metastatic melanoma* (*Nat Med*, 2019; DOI: [10.1038/s41591-019-0654-5](https://doi.org/10.1038/s41591-019-0654-5)).
* **GEO Accession**: dbGaP phs000452.v3.p1 | **cBioPortal Identifier**: `skcm_dfci_2019`
* **Cancer Type**: Metastatic Melanoma.
* **Treatment**: Anti-PD-1 monotherapy (Nivolumab or Pembrolizumab).
* **Sample Size**: **$N = 122$** pre-treatment biopsies.
* **Response Breakdown**: 48 Responders (39.3%) / 74 Non-Responders (60.7%).
* **Disk Location**: `scratch/lair/CBioPortalDataset-Liu-iAtlas`

### D. Gide et al. 2019 (`Gide-iAtlas`)
* **Publication**: *Distinct Immune Cell Populations Define Response to Anti-PD-1 Monotherapy and Anti-PD-1/Anti-CTLA-4 Combined Therapy* (*Cancer Cell*, 2019; DOI: [10.1016/j.ccell.2019.01.003](https://doi.org/10.1016/j.ccell.2019.01.003)).
* **ENA Accession**: PRJEB23709 | **cBioPortal Identifier**: `skcm_gide_2019`
* **Cancer Type**: Metastatic Melanoma.
* **Treatment**: Anti-PD-1 monotherapy (Nivolumab/Pembro) OR Combined Anti-PD-1 + Anti-CTLA-4 (Nivolumab + Ipilimumab).
* **Sample Size**: **$N = 91$** baseline biopsies.
* **Response Breakdown**: 63 Responders (69.2%) / 28 Non-Responders (30.8%).
* **Disk Location**: `scratch/lair/CBioPortalDataset-Gide-iAtlas`

### E. Rosenberg et al. 2016 / IMvigor210 (`Rosenberg-iAtlas`)
* **Publication**: *Atezolizumab in patients with locally advanced and metastatic urothelial carcinoma* (*Lancet*, 2016; DOI: [10.1016/S0140-6736(16)00561-4](https://doi.org/10.1016/S0140-6736(16)00561-4); *Nature*, 2018; DOI: [10.1038/nature25501](https://doi.org/10.1038/nature25501)).
* **cBioPortal Identifier**: `blca_imvigor210_core`
* **Cancer Type**: Locally Advanced / Metastatic Urothelial Bladder Carcinoma.
* **Treatment**: Anti-PD-L1 monotherapy (Atezolizumab).
* **Sample Size**: **$N = 347$** bulk RNA-seq samples.
* **Response Breakdown**: 68 Responders (19.6%) / 279 Non-Responders (80.4%).
* **Disk Location**: `scratch/lair/CBioPortalDataset-Rosenberg-iAtlas`

### F. Padron et al. 2022 (`Padron-iAtlas`)
* **Publication**: *Sotigalimab and/or nivolumab with chemotherapy in first-line metastatic pancreatic cancer* (*Nat Med*, 2022; DOI: [10.1038/s41591-022-01829-9](https://doi.org/10.1038/s41591-022-01829-9)).
* **Cancer Type**: Metastatic Pancreatic Ductal Adenocarcinoma (PDAC).
* **Treatment**: Anti-PD-1 (Nivolumab) + CD40 agonist (Sotigalimab) + Gemcitabine/Nab-paclitaxel.
* **Sample Size**: **$N = 93$** pre-treatment samples.
* **Response Breakdown**: 16 Responders (17.2%) / 77 Non-Responders (82.8%).
* **Disk Location**: `scratch/lair/CBioPortalDataset-Padron-iAtlas`

### G. Anders et al. 2021 (`Anders-iAtlas`)
* **Publication**: *Immune Checkpoint Blockade in Triple-Negative Breast Cancer* (*Clin Cancer Res*, 2021).
* **Cancer Type**: Metastatic Triple-Negative Breast Cancer (TNBC).
* **Treatment**: Anti-PD-L1 (Atezolizumab) + Nab-paclitaxel.
* **Sample Size**: **$N = 31$** pre-treatment biopsies.
* **Response Breakdown**: 11 Responders (35.5%) / 20 Non-Responders (64.5%).
* **Disk Location**: `scratch/lair/CBioPortalDataset-Anders-iAtlas`

### H. McDermott et al. 2018 (`McDermott-iAtlas`)
* **Publication**: *Clinical activity and molecular correlates of response to atezolizumab alone or in combination with bevacizumab in metastatic renal cell carcinoma* (*Nat Med*, 2018; DOI: [10.1038/s41591-018-0053-3](https://doi.org/10.1038/s41591-018-0053-3)).
* **cBioPortal Identifier**: `rcc_immotion150`
* **Cancer Type**: Metastatic Clear Cell Renal Cell Carcinoma (ccRCC).
* **Treatment**: Anti-PD-L1 (Atezolizumab) +/- Anti-VEGF (Bevacizumab) vs. Sunitinib.
* **Sample Size**: **$N = 263$** pre-treatment samples.
* **Response Breakdown**: 97 Responders (36.9%) / 166 Non-Responders (63.1%).
* **Disk Location**: `scratch/lair/CBioPortalDataset-McDermott-iAtlas`

### I. Choueiri et al. 2016 (`Choueiri-iAtlas`)
* **Publication**: *Immunogenomic Analysis for Biomarker Discovery in Nivolumab-treated Renal Cell Carcinoma* (*Cancer Immunol Res*, 2016; DOI: [10.1158/2326-6066.CIR-16-0048](https://doi.org/10.1158/2326-6066.CIR-16-0048)).
* **Cancer Type**: Metastatic Clear Cell Renal Cell Carcinoma (ccRCC).
* **Treatment**: Anti-PD-1 monotherapy (Nivolumab).
* **Sample Size**: **$N = 16$** pre-treatment biopsies.
* **Response Breakdown**: 6 Responders (37.5%) / 10 Non-Responders (62.5%).
* **Disk Location**: `scratch/lair/CBioPortalDataset-Choueiri-iAtlas`

---

## 3. Response Harmonization & Clinical Schema

Across all 9 cohorts, the response variable is standardized as:
```python
# scripts/sade_feldman_deconv_validation/02_deconvolute_iatlas.py:75-92
def binarize_response(resp_str: str) -> int:
    """Standardizes RECIST and clinical response annotations to binary 1/0."""
    clean = str(resp_str).strip().upper()
    if clean in ("CR", "PR", "RESPONDER", "RESPONSE", "YES", "TRUE", "1", "DCB"):
        return 1
    elif clean in ("PD", "SD", "NON-RESPONDER", "NON_RESPONDER", "NR", "NO", "FALSE", "0", "NDB"):
        return 0
    else:
        return np.nan
```

### Strata Definitions:
1. **`Melanoma` Stratum**: Combines all 4 melanoma cohorts (`Hugo`, `Riaz`, `Liu`, `Gide`) totaling **$N = 347$ melanoma patients** (159 Responders, 188 Non-Responders).
2. **`Pan-Cancer` Stratum**: Combines all 9 cohorts across Melanoma, Bladder, Pancreatic, Breast, and Renal Cell Carcinoma totaling **$N = 1{,}097$ patients** (357 Responders, 740 Non-Responders).
3. **Cohort-Specific Strata**: Evaluated individually (`Cohort_Hugo-iAtlas`, `Cohort_Riaz-iAtlas`, etc.) to calculate cross-cohort mean ROC-AUC.
