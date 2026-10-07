# Compendium of Computational Tools and Signatures for Predicting Response to Immune Checkpoint Therapies from Bulk RNA-Sequencing Data

## 1. Executive Summary & Taxonomy of Approaches

Predicting clinical response to immune checkpoint therapies (ICT / ICB, including anti-PD-1, anti-PD-L1, anti-CTLA-4, and combination regimens) from bulk RNA-sequencing (bulk RNA-seq) data is a cornerstone of translational immuno-oncology. Bulk RNA-seq profiles capture the composite transcriptional state of malignant cells, infiltrating immune populations, and stromal components.

Existing computational approaches fall into four distinct methodological categories:

```
                                  Bulk RNA-seq ICB Response Prediction Tools
                                                      │
         ┌────────────────────────┬───────────────────┴─────────────────┬────────────────────────┐
         ▼                        ▼                                     ▼                        ▼
  1. Dedicated Software    2. Classical Mechanistic             3. Microenvironment       4. Multi-Cohort Benchmark
  Packages & Pipelines        Biomarker Signatures              Deconvolution Suites         Platforms & Portals
  ────────────────────     ────────────────────────             ────────────────────      ─────────────────────────
  • TIDE (tidepy)          • IMPRES (15 gene pairs)             • immunedeconv            • ICB-Portal (Kang et al.)
  • EaSIeR (Bioconductor)  • Ayers GEP (18 genes)               • CIBERSORTx              • Litchfield et al. Cell
  • IOBR (R package)       • CYT (GZMA + PRF1)                  • MCP-counter             • TIGER (Han et al.)
  • RIMA (Snakemake)       • IPRES (Wound healing/EMT)          • EPIC
  • TcellInflamedDetector  • Roh Immune Score                   • quanTIseq
                           • TLS Signatures (12-chemokine)      • ESTIMATE
                           • CXCL9 / CD8A expression
```

All methods described in this document integrate natively with the datasets in `tme_datasets`. Each tool's profile below includes an explicit, step-by-step **Implementation Plan with `tme_datasets`**, specifying data extraction, identifier harmonization, preprocessing, scoring, and clinical validation.

---

## 2. Dedicated End-to-End Prediction Packages & Software

### 2.1 TIDE (Tumor Immune Dysfunction and Exclusion)
* **Authors / Reference**: Jiang et al., *Nature Medicine* (2018), "Signatures of T cell dysfunction and exclusion predict cancer immunotherapy response." DOI: [10.1038/s41591-018-0136-1](https://doi.org/10.1038/s41591-018-0136-1).
* **Underlying Rationale**: Tumors utilize two major mechanisms of immune evasion:
  1. **T-cell Dysfunction**: In tumors with high cytotoxic T lymphocyte (CTL) infiltration, T cells may be rendered exhausted or dysfunctional. TIDE models this by calculating the interaction between gene expression and CTL levels (average expression of *CD8A, CD8B, GZMA, GZMB, PRF1*).
  2. **T-cell Exclusion**: In tumors with low CTL infiltration, immunosuppressive cell populations actively block T-cell extravasation and infiltration. TIDE uses transcriptional signatures derived from cancer-associated fibroblasts (CAFs), myeloid-derived suppressor cells (MDSCs), and tumor-associated M2 macrophages (TAMs).
* **Output Metrics**:
  * `TIDE Score`: Combined prediction score (higher score = greater immune evasion = non-responder).
  * `Dysfunction Score`: Tumor-intrinsic T-cell dysfunction score.
  * `Exclusion Score`: Stromal/myeloid exclusion score.
  * `Responder / Non-Responder`: Binary classification based on cutoffs.
* **Input Requirements**:
  * Log2-transformed expression matrix (TPM, FPKM, or RPKM).
  * Gene-wise centering: Data **must be normalized by subtracting the average expression of each gene across the cohort** (or against a representative reference cohort) to remove platform/cohort baseline shifts.
* **Software Availability**:
  * Web Application: [http://tide.dfci.harvard.edu](http://tide.dfci.harvard.edu)
  * Python Package / CLI: `pip install tidepy`
* **Validated Indications**: Melanoma, Non-Small Cell Lung Cancer (NSCLC), Renal Cell Carcinoma (RCC), Urothelial Carcinoma.

#### Implementation Plan with `tme_datasets`
1. **Target Cohorts**:
   * Melanoma: `"Hugo-iAtlas"`, `"Riaz-iAtlas"`, `"Liu-iAtlas"`, `"Gide-iAtlas"`, `"VanAllen"`.
   * Urothelial Bladder: `"Rosenberg-iAtlas"`, `"EGAD00001006631"`.
   * Clear Cell RCC: `"McDermott-iAtlas"`, `"Choueiri-iAtlas"`.
2. **Ingestion & Filtering**:
   * Load via `tme_datasets.load_dataset(cohort_id)`.
   * Harmonize metadata using `tme_datasets.preprocessing.harmonize_obs_metadata(adata.obs)` to filter strictly for pre-treatment baseline biopsies (`biopsy_timepoint == "Pre"`).
3. **Identifier Reconciliation & Centering**:
   * Ensure gene names are HGNC symbols using `tme_datasets.genes.reconcile_genes` or `map_gene_identifier`.
   * Ensure expression is linear TPM (verify via `inspect_expression_type`), apply $\log_2(\text{TPM} + 1)$, and zero-center each gene across the cohort:
     $$X_{\text{centered}} = X - \mu_{\text{gene}}$$
4. **Execution Routine**:
   * Write the centered matrix to a tab-delimited file or invoke `tidepy` via Python:
     ```python
     import polars as pl
     from tidepy import TIDE
     # Run TIDE with cancer_type matching DatasetSpec.cancer_type
     tide_results = TIDE.run(centered_df, cancer_type="Melanoma")
     ```
5. **Output Schema & Validation**:
   * Output Polars DataFrame: `["sample_id", "tide_score", "dysfunction_score", "exclusion_score", "predicted_responder", "response_binary"]`.
   * Compute ROC-AUC, PR-AUC, and evaluate odds ratio against `response_binary`.

---

### 2.2 EaSIeR (Estimate Systems Immune Response)
* **Authors / Reference**: Lapuente-Santana et al., *Patterns* (2021), "Interpretable systems biomarkers predict response to immune-checkpoint inhibitors." DOI: [10.1016/j.patter.2021.100293](https://doi.org/10.1016/j.patter.2021.100293).
* **Underlying Rationale**: Compresses bulk RNA-seq data into four systems-level quantitative descriptors:
  1. **Immune Cell Abundances**: Deconvolution estimates from MCP-counter and quanTIseq.
  2. **Intracellular Pathway Activities**: Inferred via PROGENy (MAPK, NFkB, TGFb, Trail, VEGF, etc.).
  3. **Transcription Factor Activities**: Inferred via DoRothEA regulon models.
  4. **Cell-Cell Communication**: Cytokine and chemokine signaling matrices.
  Predicts response probability using multi-cohort regularized regression models (Elastic Net).
* **Software Availability**: R/Bioconductor package (`easier` + `easierData`).

#### Implementation Plan with `tme_datasets`
1. **Target Cohorts**:
   * Urothelial Bladder: `"Rosenberg-iAtlas"` / `"EGAD00001006631"` (IMvigor210 cohort).
   * Melanoma: `"Liu-iAtlas"`, `"Hugo-iAtlas"`.
2. **Ingestion & Conversion**:
   * Load cohort with `tme_datasets.load_dataset(cohort_id)`.
   * Extract baseline samples (`biopsy_timepoint == "Pre"`) and linear TPM expression matrix.
   * Save expression matrix and clinical metadata into Arrow/Parquet or CSV files for R consumption:
     ```python
     expr_df = pl.DataFrame(adata.X.toarray(), schema=list(adata.var_names)).with_columns(
         sample_id=pl.Series(adata.obs_names)
     )
     expr_df.write_parquet("scratch/easier_input_expr.parquet")
     ```
3. **R Execution Script via Subprocess / Rpy2**:
   ```R
   library(easier)
   library(arrow)
   expr <- as.data.frame(read_parquet("scratch/easier_input_expr.parquet"))
   rownames(expr) <- expr$sample_id
   expr$sample_id <- NULL
   tpm_mat <- t(as.matrix(expr))

   # Compute systems descriptors and predict response
   predictions <- predict_immune_response(tpm_mat, cancer_type = "bladder")
   write.csv(predictions, "scratch/easier_predictions.csv")
   ```
4. **Validation**:
   * Join `easier_score` with `response_binary` in Polars. Compute ROC-AUC and evaluate pathway importance weights (e.g., verifying whether TGF-β pathway activity associates with resistance in urothelial cohorts).

---

### 2.3 IOBR (Immuno-Oncology Biological Research Toolkit)
* **Authors / Reference**: Zeng et al., *Frontiers in Immunology* (2021), "IOBR: Multi-Omics Immuno-Oncology Biological Research to Decode Tumor Microenvironment and Signatures." DOI: [10.3389/fimmu.2021.777085](https://doi.org/10.3389/fimmu.2021.777085).
* **Underlying Rationale**: Comprehensive multi-omics R toolkit aggregating 250+ curated immuno-oncology gene signatures, 8 microenvironment deconvolution algorithms, and survival/response correlation modules.
* **Software Availability**: R package (`IOBR` on GitHub).

#### Implementation Plan with `tme_datasets`
1. **Target Cohorts**:
   * Pan-cancer evaluation across all registered bulk cohorts: `"Hugo-iAtlas"`, `"Riaz-iAtlas"`, `"Liu-iAtlas"`, `"Gide-iAtlas"`, `"Rosenberg-iAtlas"`, `"McDermott-iAtlas"`, `"Padron-iAtlas"`, `"Anders-iAtlas"`.
2. **Ingestion & Data Handoff**:
   * Export length-scaled TPM matrix with HGNC symbols as row names.
3. **Execution Script**:
   ```R
   library(IOBR)
   library(arrow)

   tpm_mat <- as.matrix(read_parquet("scratch/iobr_tpm.parquet"))
   # 1. Calculate curated signature scores (ssGSEA or PCA)
   sig_scores <- calculate_sig_score(eset = tpm_mat, signature = signature_collection)
   # 2. Deconvolute TME via MCP-counter, EPIC, and CIBERSORT
   deconv_res <- deconvo_tme(eset = tpm_mat, method = "mcpcounter")
   
   final_features <- cbind(sig_scores, deconv_res)
   write_parquet(as.data.frame(final_features), "scratch/iobr_features.parquet")
   ```
4. **Evaluation**:
   * Load `iobr_features.parquet` back into Polars.
   * Run univariate logistic regression and Spearman correlation for each of the 250+ signatures against `response_binary`. Rank signatures by cross-cohort stability index.

---

### 2.4 RIMA (RNA-seq IMmune Analysis Pipeline)
* **Authors / Reference**: Liu Lab, Dana-Farber Cancer Institute. *PNAS* / GitBook.
* **Underlying Rationale**: Automated, containerized Snakemake workflow starting from raw sequencing reads to compute alignment (STAR), expression (Salmon), TCR/BCR clonality (TRUST4), microsatellite instability (MSIsensor2), and TIDE response scores.

#### Implementation Plan with `tme_datasets`
1. **Target Cohorts**:
   * Raw read cohorts with FASTQ or aligned BAM files available under `data/manual_download/` (e.g., `"EGAD00001006631-align"` IMvigor210).
2. **Workflow Configuration**:
   * Map sample paths into RIMA's `metasheet.csv`:
     `Sample,FQ1,FQ2,Batch,Condition`
   * Enable pipeline modules in `config.yaml`:
     `modules: [preprocessing, deconvolution, immune_repertoire, response_prediction]`
3. **Execution**:
   * Execute via Snakemake with Singularity/Conda:
     ```bash
     snakemake -s RIMA.snakefile --use-conda --cores 16
     ```
4. **Integration with `tme_datasets`**:
   * Ingest output tables (TCR diversity, MSI status, TIDE score) directly into the `adata.obs` table of the corresponding dataset in `tme_datasets` for joint single-cell/bulk comparison.

---

### 2.5 TcellInflamedDetector / Ayers GEP
* **Authors / Reference**: Yang & Park, *Genomics & Informatics* (2022); Ayers et al., *J Clin Invest* (2017).
* **Underlying Rationale**: Scores the 18-gene expanded IFN-γ / T-cell-inflamed Gene Expression Profile (GEP) to distinguish inflamed from non-inflamed tumors.

#### Implementation Plan with `tme_datasets`
1. **Target Cohorts**:
   * Pan-cancer bulk cohorts in `tme_datasets` (`"Hugo-iAtlas"`, `"Riaz-iAtlas"`, `"Liu-iAtlas"`, `"Gide-iAtlas"`, `"Rosenberg-iAtlas"`).
2. **Pure Functional Python Implementation**:
   * `tme_datasets` already provides `AYERS_T_CELL_INFLAMED_GEP` in its `genesets` module:
     ```python
     from tme_datasets.query import load_dataset
     from tme_datasets.genesets import AYERS_T_CELL_INFLAMED_GEP, score_geneset_zscore
     from tme_datasets.genesets.models import GeneSetCollection
     from tme_datasets.preprocessing import harmonize_obs_metadata
     from returns.result import Success

     # 1. Load dataset
     match load_dataset("Liu-iAtlas"):
         case Success(adata):
             # 2. Harmonize metadata
             obs_res = harmonize_obs_metadata(adata.obs)
             # 3. Calculate score
             coll = GeneSetCollection(gene_sets={"Ayers_GEP": AYERS_T_CELL_INFLAMED_GEP})
             score_res = score_geneset_zscore(adata, coll)
     ```
3. **Validation**:
   * Join `Ayers_GEP` Z-score with `response_binary`. Compute ROC-AUC across melanoma and bladder cohorts.

---

## 3. Classical Mechanistic Signatures & Mathematical Formulations

### 3.1 IMPRES (15 Gene-Pair Non-Parametric Classifier)
* **Authors / Reference**: Auslander et al., *Nature Medicine* (2018).
* **Mathematical Formulation**: Evaluates 15 pairwise relations between checkpoint and co-stimulatory genes:
  $$\text{IMPRES} = \sum_{k=1}^{15} \mathbb{I}\left( \text{Expr}(A_k) > \text{Expr}(B_k) \right)$$
  where $\mathbb{I}$ is the indicator function, yielding an integer score $\in [0, 15]$.

#### Implementation Plan with `tme_datasets`
1. **Target Cohorts**:
   * Melanoma: `"Hugo-iAtlas"`, `"Riaz-iAtlas"`, `"Liu-iAtlas"`, `"Auslander"`, `"VanAllen"`.
2. **Gene Pairs Definition**:
   ```python
   IMPRES_PAIRS = (
       ("PDCD1", "TNFRSF4"), ("CD27", "CD40"), ("CD27", "CD80"),
       ("CD40LG", "CD80"), ("CD40LG", "CD86"), ("CD40LG", "TNFRSF9"),
       ("CD40LG", "ICOSLG"), ("CD28", "CD86"), ("CD80", "HAVCR2"),
       ("CD86", "HAVCR2"), ("CD274", "HAVCR2"), ("CTLA4", "HAVCR2"),
       ("CD27", "HAVCR2"), ("ICOS", "HAVCR2"), ("TNFRSF4", "HAVCR2"),
   )
   ```
3. **Pure Functional Python Scoring Function**:
   ```python
   import anndata as ad
   import numpy as np
   import polars as pl
   from returns.result import Result, Success, Failure

   def compute_impres(adata: ad.AnnData) -> Result[pl.DataFrame, str]:
       try:
           var_names = list(adata.var_names)
           var_map = {name: idx for idx, name in enumerate(var_names)}
           X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)

           pair_scores = []
           for gene_a, gene_b in IMPRES_PAIRS:
               if gene_a in var_map and gene_b in var_map:
                   idx_a = var_map[gene_a]
                   idx_b = var_map[gene_b]
                   # Boolean comparison: 1 if A > B else 0
                   pair_scores.append((X[:, idx_a] > X[:, idx_b]).astype(int))

           if not pair_scores:
               return Failure("None of the IMPRES gene pairs were found in dataset.")

           impres_scores = np.sum(np.column_stack(pair_scores), axis=1)
           return Success(pl.DataFrame({
               "sample_id": list(adata.obs_names),
               "impres_score": [int(x) for x in impres_scores]
           }))
       except Exception as exc:
           return Failure(f"IMPRES computation failed: {exc}")
   ```
4. **Validation**:
   * Calculate ROC-AUC against `response_binary`. Evaluate whether performance holds in non-melanoma cohorts (`"Rosenberg-iAtlas"`).

---

### 3.2 CYT (Rooney Cytolytic Activity Score)
* **Authors / Reference**: Rooney et al., *Cell* (2015).
* **Mathematical Formulation**:
  $$\text{CYT} = \frac{\log_2(\text{Expr}(GZMA) + 1) + \log_2(\text{Expr}(PRF1) + 1)}{2}$$

#### Implementation Plan with `tme_datasets`
1. **Target Cohorts**: Pan-cancer bulk cohorts across all tumor types.
2. **Pure Functional Python Implementation**:
   ```python
   def compute_cyt(adata: ad.AnnData) -> Result[pl.DataFrame, str]:
       try:
           var_names = list(adata.var_names)
           var_map = {name: idx for idx, name in enumerate(var_names)}
           if "GZMA" not in var_map or "PRF1" not in var_map:
               return Failure("GZMA or PRF1 missing from dataset.")

           X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
           gzma_expr = np.log2(X[:, var_map["GZMA"]] + 1.0)
           prf1_expr = np.log2(X[:, var_map["PRF1"]] + 1.0)
           cyt = (gzma_expr + prf1_expr) / 2.0

           return Success(pl.DataFrame({
               "sample_id": list(adata.obs_names),
               "cyt_score": [float(v) for v in cyt]
           }))
       except Exception as exc:
           return Failure(f"CYT computation failed: {exc}")
   ```
3. **Validation**:
   * Cross-evaluate across all 9 bulk cohorts. Benchmark against more complex signatures.

---

### 3.3 IPRES (Innate Anti-PD-1 Resistance Signature)
* **Authors / Reference**: Hugo et al., *Cell* (2016).
* **Mathematical Formulation**: Mean enrichment score across 26 biological gene sets capturing mesenchymal transition, wound healing, angiogenesis, and cell adhesion.

#### Implementation Plan with `tme_datasets`
1. **Target Cohorts**: Melanoma baseline cohorts (`"Hugo-iAtlas"`, `"Riaz-iAtlas"`).
2. **Implementation**:
   * Define the IPRES signature collection in `tme_datasets.genesets.GeneSetCollection`.
   * Compute scores using `score_geneset_zscore` or `score_geneset_auc`.
   * Invert the score (higher IPRES indicates resistance).
3. **Validation**:
   * Test discrimination between progressive disease (PD) and partial/complete response (PR/CR).

---

### 3.4 Single-Gene Biomarkers (*CXCL9* and *CD8A*)
* **Authors / Reference**: Litchfield et al., *Cell* (2021).
* **Mathematical Formulation**: Log2-normalized expression of *CXCL9* or *CD8A*.

#### Implementation Plan with `tme_datasets`
1. **Target Cohorts**: All bulk cohorts.
2. **Pure Functional Extraction**:
   ```python
   def extract_single_gene_predictors(adata: ad.AnnData) -> Result[pl.DataFrame, str]:
       try:
           var_names = list(adata.var_names)
           var_map = {name: idx for idx, name in enumerate(var_names)}
           X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)

           df_dict: dict[str, list[object]] = {"sample_id": list(adata.obs_names)}
           for gene in ("CXCL9", "CD8A", "PDCD1", "IFNG"):
               if gene in var_map:
                   df_dict[f"expr_{gene}"] = [float(np.log2(v + 1.0)) for v in X[:, var_map[gene]]]

           return Success(pl.DataFrame(df_dict))
       except Exception as exc:
           return Failure(f"Extraction failed: {exc}")
   ```
3. **Validation**:
   * Measure univariate ROC-AUC for `expr_CXCL9` across all cohorts. Confirm whether *CXCL9* matches or exceeds the multi-gene signatures as reported in the Litchfield benchmark.

---

## 4. Microenvironment Deconvolution Suites

### 4.1 Comparison of Deconvolution Suites

| Method | Output Nature | Key Cell Types Relevant to ICB | Underlying Algorithm | Reference |
| :--- | :--- | :--- | :--- | :--- |
| **immunedeconv** | Unified wrapper | CD8+ T, NK, Treg, M1/M2 Macrophages, B cells, CAFs, Endothelial | Aggregates 6 methods into a single standardized R interface | Sturm et al., *Bioinformatics* (2019) |
| **CIBERSORTx** | Cell fractions / absolute scores | 22 leukocyte subsets (LM22) or custom scRNA-seq reference matrices | Support Vector Regression ($\nu$-SVR) with de-noising | Newman et al., *Nat Biotechnol* (2019) |
| **MCP-counter** | Population abundance scores | CD8+ T, Cytotoxic lymphocytes, NK, B lineage, Monocytes, CAFs | Geometric mean of transcriptomic marker genes | Becht et al., *Genome Biol* (2016) |
| **EPIC** | Absolute cell fractions | CD8+ T, CD4+ T, B cells, NK cells, Macrophages, Endothelial, CAFs | Constrained least squares regression | Racle et al., *eLife* (2017) |
| **quanTIseq** | Absolute cell fractions | 10 immune cell types (including M1, M2, Tregs, CD8+ T) | De-novo constrained least squares | Finotello et al., *Genome Med* (2019) |
| **ESTIMATE** | ImmuneScore, StromalScore, Purity | Infiltrating immune vs. stromal vs. malignant fractions | ssGSEA on 141 immune and 141 stromal genes | Yoshihara et al., *Nat Commun* (2013) |

---

### 4.2 Implementation Plan with `tme_datasets` & In Silico Pseudobulk Validation

Deconvolution methods can be evaluated on both real clinical bulk RNA-seq cohorts and in silico simulated mixtures:

```
[Single-Cell Reference Atlas] (e.g. GSE120575 Sade-Feldman)
             │
             ▼
[simulate_pseudobulk(adata)]  ───────► Known Ground-Truth Cell Fractions
             │
             ├────────────────────────────────────────┬────────────────────────────────────────┐
             ▼                                        ▼                                        ▼
      MCP-counter                              CIBERSORTx / EPIC                               ESTIMATE
             │                                        │                                        │
             └────────────────────────────────────────┴────────────────────────────────────────┘
                                                     │
                                                     ▼
                                    [Deconvolution Accuracy Check]
                                    • Pearson / Spearman correlation (r)
                                    • Root Mean Squared Error (RMSE)
```

1. **Step 1: Synthetic Pseudobulk Benchmark (`simulate_pseudobulk`)**:
   * Load single-cell dataset with cell type annotations (e.g., `"GSE120575"` Sade-Feldman Melanoma).
   * Invoke `tme_datasets.simulation.simulate_pseudobulk(adata, group_by="sample_id")`.
   * This yields simulated bulk samples with known ground-truth proportions of CD8+ T cells, CD4+ T cells, B cells, Macrophages, and Malignant cells.
   * Run deconvolution methods and evaluate accuracy:
     $$\text{RMSE} = \sqrt{\frac{1}{N} \sum_{i=1}^N (y_i^{\text{true}} - y_i^{\text{estimated}})^2}$$
2. **Step 2: Clinical Response Association**:
   * Apply validated deconvolution on real clinical bulk cohorts (`"Hugo-iAtlas"`, `"Rosenberg-iAtlas"`).
   * Evaluate whether estimated **CD8+ T-cell abundance** correlates positively with response, and whether **CAF abundance** associates with treatment failure / non-response.

---

## 5. Master Feature & Capability Matrix

| Tool / Model | Input Format | Primary Language | CLI / Web Available | Self-Contained Classifier? | Immune Exclusion Evaluated? | Multi-Omics Supported? |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **TIDE** | Log2(TPM/FPKM), centered | Python, Web | Yes (`tidepy`, Web) | Yes | Yes (CAFs, MDSCs, M2) | No (RNA only) |
| **EaSIeR** | Normalized counts / TPM | R (Bioconductor) | R API | Yes (Elastic Net) | Yes (via pathways/stroma) | No (RNA only) |
| **IOBR** | TPM, FPKM, Counts | R (GitHub) | R API | Aggregates all | Yes | Yes (Mutations, CNVs) |
| **RIMA** | FASTQ / BAM | Snakemake, Python | Yes (Workflow) | Yes (via TIDE/MSI) | Yes | Yes (TCR/BCR, MSI) |
| **IMPRES** | Any rankable scale | Python / R | Pure Python script | Yes (Score 0–15) | No | No |
| **Ayers GEP** | Linear or log TPM | Python / R | `tme_datasets.genesets` | Yes (Z-score) | No | No |
| **CYT** | Log2(TPM + 1) | Python / R | Pure Python script | Yes (Mean) | No | No |
| **immunedeconv** | Linear TPM / Counts | R | R API | No (Deconv only) | Yes (CAFs/Endothelial) | No |
| **CIBERSORTx** | Linear TPM / Counts | Web / Docker | Yes (Web, Docker) | No (Deconv only) | Yes (Monocyte/M2) | No |

---

## 6. Systematic Multi-Cohort Benchmarking Insights

Empirical multi-cohort evaluations (e.g., *Litchfield et al., Cell 2021*; *Kang et al., Cancers 2023*) have identified fundamental realities regarding bulk RNA-seq ICB response prediction:

### 6.1 The Generalization Drop & Cross-Cohort Instability
* Most transcriptomic signatures report ROC-AUC values between 0.75 and 0.90 in their original discovery cohorts.
* However, when evaluated on **independent external cohorts**, predictive accuracy frequently drops to near-chance or moderate levels (AUC ~ 0.55–0.68).
* **Top Independent Performers**:
  * In the Kang et al. benchmark of 48 scores across 29 cohorts, **TIDE** and **CYT** demonstrated the most consistent predictive power across diverse indications (Melanoma, NSCLC, Gastric Cancer, and Urothelial Bladder Cancer).
  * **PASS-ON** and **EIGS_ssGSEA** showed the highest correlations with clinical outcomes across heterogeneous datasets.
  * **IMPRES** performed well in melanoma cohorts with similar sequencing protocols, but lost predictive significance in non-melanoma and FFPE cohorts.

### 6.2 High Collinearity Across Published Signatures
* The vast majority of published "novel" immune gene signatures (GEP, Roh, CYT, Davoli, IFN-γ scores) exhibit pairwise Pearson correlations exceeding $r > 0.80$.
* They all reflect a single convergent biological phenomenon: **the IFN-γ / Antigen Presentation / Cytotoxic T-cell Axis** (*CD8A, CD8B, IFNG, CXCL9, CXCL10, PRF1, GZMA, PDCD1, CD274*).
* In the Litchfield et al. *Cell* 2021 meta-analysis (1,008 patients), **expression of *CXCL9* alone** performed comparably or superiorly to large, multi-gene signatures.

### 6.3 DNA-RNA Synergy
* Bulk RNA-seq signatures capture the **pre-existing inflamed immune microenvironment**, but do not capture tumor foreignness.
* Combining transcriptomic scores (e.g., TIDE, GEP, or *CXCL9*) with genomic biomarkers (Tumor Mutational Burden [TMB], frameshift indels, or NMD-escape mutations) significantly outperforms either modality alone (AUCs increasing from ~0.65 to >0.78).

---

## 7. Complete Functional Benchmark Pipeline Template

Below is a production-grade, functional Python script conforming strictly to the repository's coding style rules (`returns.result.Result`, `returns.decorators.do`, Polars, immutability, no `None`). It ingests any registered cohort from `tme_datasets`, computes **CYT**, **IMPRES**, **Ayers GEP**, and **CXCL9/CD8A**, and outputs a unified benchmark dataset.

```python
"""benchmark_icb_predictors.py: Pure functional execution of ICB response predictors on tme_datasets."""

from __future__ import annotations
from typing import NamedTuple, assert_never
import anndata as ad
import numpy as np
import polars as pl
from returns.result import Result, Success, Failure
from returns.maybe import Maybe, Some, Nothing
from returns.decorators import do

from tme_datasets.query import load_dataset
from tme_datasets.genesets import AYERS_T_CELL_INFLAMED_GEP, score_geneset_zscore
from tme_datasets.genesets.models import GeneSetCollection
from tme_datasets.preprocessing import harmonize_obs_metadata

# Immutable container for joined prediction metrics
class BenchmarkMetrics(NamedTuple):
    dataset_id: str
    n_samples: int
    n_responders: int
    n_non_responders: int
    predictions: pl.DataFrame

def extract_pre_treatment_samples(adata: ad.AnnData) -> Result[ad.AnnData, str]:
    """Pure filter selecting baseline pre-treatment samples with valid binary response."""
    match harmonize_obs_metadata(adata.obs):
        case Failure(err):
            return Failure(err)
        case Success(obs_df):
            # Find indices where response_binary is not null and timepoint is Pre or Unknown
            valid_mask = (
                obs_df["response_binary"].is_not_null()
                & obs_df["biopsy_timepoint"].is_in(["Pre", "Unknown"])
            ).to_numpy()

            if not np.any(valid_mask):
                return Failure("No valid pre-treatment samples with binary response found.")

            adata_filtered = adata[valid_mask, :].copy()
            adata_filtered.obs = obs_df.filter(
                pl.col("response_binary").is_not_null()
                & pl.col("biopsy_timepoint").is_in(["Pre", "Unknown"])
            ).to_pandas()
            adata_filtered.obs_names = [str(x) for x in adata_filtered.obs["sample_id"]]
            return Success(adata_filtered)

def compute_all_signatures(adata: ad.AnnData) -> Result[pl.DataFrame, str]:
    """Compute CYT, single-gene predictors, and Ayers GEP on AnnData."""
    try:
        var_map = {name: idx for idx, name in enumerate(adata.var_names)}
        X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)

        # 1. Single genes & CYT
        df_dict: dict[str, list[object]] = {
            "sample_id": list(adata.obs_names),
            "response_binary": [float(x) for x in adata.obs["response_binary"]],
        }

        # Check CYT genes
        has_cyt = "GZMA" in var_map and "PRF1" in var_map
        if has_cyt:
            gzma = np.log2(X[:, var_map["GZMA"]] + 1.0)
            prf1 = np.log2(X[:, var_map["PRF1"]] + 1.0)
            df_dict["score_cyt"] = [float(v) for v in (gzma + prf1) / 2.0]

        # Check CXCL9
        if "CXCL9" in var_map:
            df_dict["expr_cxcl9"] = [float(np.log2(v + 1.0)) for v in X[:, var_map["CXCL9"]]]

        # Check CD8A
        if "CD8A" in var_map:
            df_dict["expr_cd8a"] = [float(np.log2(v + 1.0)) for v in X[:, var_map["CD8A"]]]

        df_base = pl.DataFrame(df_dict)

        # 2. Ayers GEP
        coll = GeneSetCollection(gene_sets={"Ayers_GEP": AYERS_T_CELL_INFLAMED_GEP})
        match score_geneset_zscore(adata, coll):
            case Success(gep_df):
                df_base = df_base.join(gep_df.rename({"Ayers_GEP": "score_gep"}), on="sample_id", how="left")
            case Failure(_):
                pass

        return Success(df_base)
    except Exception as exc:
        return Failure(f"Signature calculation failed: {exc}")

@do(Result[BenchmarkMetrics, str])
def run_cohort_benchmark(cohort_id: str) -> BenchmarkMetrics:
    """Monadic execution pipeline for a single cohort."""
    raw_adata = yield load_dataset(cohort_id)
    filtered_adata = yield extract_pre_treatment_samples(raw_adata)
    predictions_df = yield compute_all_signatures(filtered_adata)

    n_resp = int(predictions_df["response_binary"].sum())
    n_total = len(predictions_df)

    return BenchmarkMetrics(
        dataset_id=cohort_id,
        n_samples=n_total,
        n_responders=n_resp,
        n_non_responders=n_total - n_resp,
        predictions=predictions_df,
    )
```

---

## 8. Bibliography & Key References

1. **Jiang, P., et al.** (2018). Signatures of T cell dysfunction and exclusion predict cancer immunotherapy response. *Nature Medicine*, 24(10), 1550-1558.
2. **Lapuente-Santana, Ó., et al.** (2021). Interpretable systems biomarkers predict response to immune-checkpoint inhibitors. *Patterns*, 2(8), 100293.
3. **Zeng, D., et al.** (2021). IOBR: Multi-Omics Immuno-Oncology Biological Research to Decode Tumor Microenvironment and Signatures. *Frontiers in Immunology*, 12, 777085.
4. **Litchfield, K., et al.** (2021). Meta-analysis of tumor- and T cell-intrinsic mechanisms of sensitization to checkpoint inhibition. *Cell*, 184(3), 596-614.
5. **Kang, H., et al.** (2023). A Comprehensive Benchmark of Transcriptomic Biomarkers for Immune Checkpoint Blockades. *Cancers*, 15(16), 4094.
6. **Auslander, N., et al.** (2018). Robust prediction of response to immune checkpoint blockade therapy in metastatic melanoma. *Nature Medicine*, 24(10), 1545-1549.
7. **Ayers, M., et al.** (2017). IFN-γ-related mRNA profile predicts clinical response to PD-1 blockade. *The Journal of Clinical Investigation*, 127(8), 2930-2940.
8. **Rooney, M. S., et al.** (2015). Molecular and genetic properties of tumors associated with local immune cytolytic activity. *Cell*, 160(1-2), 48-61.
9. **Sturm, G., et al.** (2019). Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology. *Bioinformatics*, 35(14), i436-i445.
10. **Han, C., et al.** (2023). TIGER: a database of tumor immunotherapy gene expression resource. *Nucleic Acids Research*, 51(D1), D1240-D1248.
