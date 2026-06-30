# Report: Mitigation of PCA Cohort-Splitting Batch Effects between TCGA and iAtlas Datasets

This report documents the diagnostic findings, root cause, and mitigation steps taken to resolve the batch effect that caused the combined TCGA and iAtlas expression data to split cleanly along the first principal component (PC1).

---

## 1. Problem Statement
In the initial combined PCA analysis of the TCGA and iAtlas datasets (using both DE genes and Bagaev immune signatures), **PC1 separated the iAtlas cohorts from the TCGA cohort** with near-perfect accuracy (Silhouette Score > 0.90).

---

## 2. Root Cause Analysis & Diagnostic Findings

### 2.1 PC1 Loading Contributions
An investigation into the gene contributions (loadings) for PC1 revealed a systematic split:
- **Positive Loadings (High in iAtlas)**: Lowly-expressed immune marker genes (e.g. `CD28`, `ICOS`, `CD80`, `ITK`) had mean values 10-100x higher in iAtlas than in TCGA (e.g., `CD28` mean was $3.5 \times 10^{-5}$ in iAtlas vs. $2.58 \times 10^{-6}$ in TCGA).
- **Negative Loadings (High in TCGA)**: Highly-expressed housekeeping/structural genes (e.g. `B2M`, `HLA-A`, `HLA-B`, `HLA-C`, `COL1A1`) had mean values 15-30x higher in TCGA than in iAtlas (e.g., `B2M` mean was $2.74 \times 10^{-3}$ in TCGA vs. $9.4 \times 10^{-5}$ in iAtlas).

### 2.2 Log-Space Normalization Discrepancy
- The raw TCGA counts were correctly processed in **linear space** before normalization.
- However, the iAtlas datasets downloaded from cBioPortal are pre-stored as **log2-transformed values** (e.g. $\log_2(\text{TPM} + 1)$). The sum of these log-scale values over a sample ranged between ~130,000 and ~220,000 (instead of $1,000,000$ for linear TPM).
- The pipeline was performing sum-to-one normalization **directly on the log2 values** (i.e. dividing the log values by the sum of log values). 
- Because division in log-space is mathematically incorrect, it severely compressed the dynamic range:
  - The ratio of highly expressed `B2M` to lowly expressed `CD28` was compressed to **2.3** in iAtlas (down from a biological ratio of **1062** in TCGA).
  - This flattening of the dynamic range made highly-expressed genes appear vastly under-expressed in iAtlas, and lowly-expressed genes appear vastly over-expressed, creating a major batch effect.

---

## 3. Mitigation Implementation

We updated the preprocessing and loading sequence in [TCGA-background.py](file:///Users/halu/Code/tme_analysis/scripts/TCGA-background.py) to process both datasets in the exact same mathematical scale:

1. **iAtlas Un-logging**:
   Introduced a custom loader `load_and_transform_iatlas_properly` which detects the file type and applies $2^y - 1$ to convert the log2-transformed values back into linear abundance values.
2. **Gene Length (RPK) Correction**:
   If the raw iAtlas file contains counts, we perform gene length correction using exon union lengths retrieved from SQLite queries against the Ensembl release 111 database.
3. **Linear Sum-to-One Normalization**:
   Linear expression values for both datasets are normalized so that their column sums equal exactly `1.0` before applying the log-transformation `log1p`. Both matrices are now in the exact same physical units:
   $$\text{log1p}\left( \frac{\text{Linear\_abundance}_i}{\sum_j \text{Linear\_abundance}_j} \right)$$

---

## 4. Results & Verification
- **Cohort Mixing**: The artificial split along PC1 was successfully mitigated. Cohorts now cluster by biological similarity (e.g., melanoma iAtlas cohorts overlap SKCM TCGA samples).
- **PCA Silhouette Score Drop**: The best K-means Silhouette Score in the PCA space fell from **0.93** (representing batch separation) to a healthy **0.55** (representing integrated cohorts).
- **File Outputs**: All updated plots and metrics reports have been successfully generated and copied to the local workspace and artifacts directory.
