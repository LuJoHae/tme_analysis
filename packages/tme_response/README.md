# `tme_response`: Multi-Omic Immune Checkpoint Therapy Response Prediction Framework

`tme_response` is a functional-first Python library designed to predict and benchmark patient response to immune checkpoint therapies (ICT / ICB, including anti-PD-1, anti-PD-L1, anti-CTLA-4) by integrating bulk transcriptomic and genomic data.

It directly depends on the workspace package `tme_datasets` for access to standardized bulk RNA-seq cohorts, single-cell reference atlases (for pseudobulk simulation), curated gene sets, and clinical annotations.

---

## 1. Key Features

* **Multi-Omics DNA + RNA Integration**:
  * Seamlessly pairs bulk transcriptomic states with genomic determinants (TMB, frameshift indels, aneuploidy, and loss-of-heterozygosity [LOH] / inactivating driver mutations like *B2M*, *JAK1*, *JAK2*, *STK11*, *KEAP1*).
* **Comprehensive Signature Suite**:
  * **CYT** (Rooney Cytolytic Activity: $\sqrt{\text{GZMA} \cdot \text{PRF1}}$).
  * **IMPRES** (Auslander 15-pair non-parametric score).
  * **Ayers Expanded GEP** (18-gene T-cell inflamed profile via `tme_datasets.genesets.AYERS_T_CELL_INFLAMED_GEP`).
  * **Single-Gene Biomarkers** (Litchfield *CXCL9*, *CD8A*, *PDCD1*, *IFNG*).
  * **IPRES** (Hugo 26-gene set innate resistance).
* **Systems Biology & Evasion Models**:
  * **TIDE** (Tumor Immune Dysfunction and Exclusion with cohort-level zero-centering).
  * **EaSIeR** and **IOBR** interop bridges.
* **Cellular Infiltration & Deconvolution**:
  * MCP-counter CD8+ T-cell, cytotoxic lymphocyte, and CAF infiltrate estimators.
  * In silico pseudobulk benchmark validation using `tme_datasets.simulation.simulate_pseudobulk`.
* **Synergy Models**:
  * Linear composite models: $S_{\text{composite}} = z(\text{RNA}) + \gamma \cdot z(\log_{10}(\text{TMB} + 1))$.
  * Biological gating logic (e.g. *CXCL9*-high + *B2M*-WT).
* **Declarative Visualization**:
  * ROC curves, PR curves, and cross-cohort benchmark comparisons in Altair, exported directly to publication-ready vector SVGs.
* **Strict Functional Programming Standards**:
  * Fully typed, pure functions, monadic error handling via `returns` (`Result` and `Maybe`), Polars dataframes, and property-based testing with `hypothesis`.

---

## 2. Quick Start Example

```python
from returns.result import Success
from tme_response import (
    load_multi_omic_cohort,
    compute_cyt_score,
    compute_impres_score,
    compute_single_gene_score,
    compute_dna_rna_composite,
    evaluate_prediction,
)

# 1. Ingest synchronized multi-omic cohort from tme_datasets
match load_multi_omic_cohort("Liu-iAtlas"):
    case Success(cohort):
        print(f"Loaded {cohort.cohort_id} with {len(cohort.sample_ids)} baseline samples.")

        # 2. Compute transcriptomic signatures
        cyt_pred = compute_cyt_score(cohort).unwrap()
        cxcl9_pred = compute_single_gene_score(cohort, "CXCL9").unwrap()

        # 3. Compute DNA-RNA synergy composite (if TMB available)
        composite_pred = compute_dna_rna_composite(cohort, cxcl9_pred).unwrap()

        # 4. Evaluate ROC-AUC and PR-AUC against RECIST binary response
        eval_res = evaluate_prediction(
            composite_pred,
            cohort.clinical_annotations,
            cohort.cohort_id,
            cohort.cancer_type,
        ).unwrap()

        print(f"Composite Model ROC-AUC: {eval_res.roc_auc:.3f}, PR-AUC: {eval_res.pr_auc:.3f}")
```

---

## 3. Architecture

```
tme_response/
├── data/           # MultiOmicCohort adapter, TMB/MAF parser, CNA/aneuploidy extractor
├── signatures/     # Pure functional CYT, IMPRES, Ayers GEP, CXCL9, IPRES
├── models/         # TIDE (zero-centering, dysfunction, exclusion), EaSIeR bridge
├── synergy/        # Composite DNA-RNA scores & mutation-gating logic
├── deconvolution/  # Infiltrate estimators & pseudobulk validation harness
├── evaluation/     # Metrics (ROC-AUC, PR-AUC, Brier score) & multi-cohort benchmark runner
└── visualization/  # Altair vector graphics with SVG export
```
