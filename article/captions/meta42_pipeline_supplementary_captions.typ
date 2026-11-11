// =============================================================================
// Importable Captions for Pipeline Parameter Optimization & Consensus-Meta42
// Supplementary Figures (S57-S60)
// =============================================================================

#let caption_supp_pipeline_envelope = [
  *Dependence of maximum achievable bulk ICB response classification accuracy on single-cell reference construction parameters across 312 empirical configurations.*
  *(a)* Empirical upper-bound performance envelope showing maximum, 90th-percentile ($Q_{90}$), 75th-percentile ($Q_{75}$), and median cross-validated response classification ROC-AUC as a function of total reference cell states ($K$, log scale), highlighting the `Consensus-Meta42` hierarchical optimum ($K = 42$, $"AUC" = 0.746$) relative to the unharmonized 9-cohort ceiling.
  *(b)* Bivariate pipeline parameter optimization surface mapping mean cohort multivariate ROC-AUC across the $9 times 8$ grid of combined single-cell datasets ($1 dots 9$) and Leiden clustering resolutions ($"res" in [0.25, 2.00]$).
  *(c)* Random Forest variance decomposition ($R^2$ permutation importance) quantifying the relative contribution of reference cell-state granularity ($54.8%$), signature matrix collinearity ($29.9%$), clustering resolution ($6.8%$), multi-cohort dataset depth ($4.7%$), and total single-cell count ($3.7%$) to overall predictive performance.
  *(d)* Stepwise hyperparameter optimization waterfall tracing the cumulative ROC-AUC progression from an unoptimized single-cohort baseline ($"res" = 0.25$, $"AUC" = 0.612$) through multi-cohort expansion ($+0.036$), high-resolution sub-clustering ($+0.056$), MiloR differential abundance state selection ($+0.025$), and Ward hierarchical consensus deduplication ($+0.017$) to reach the `Consensus-Meta42` peak ($"AUC" = 0.746$).
]

#let caption_supp_parameter_mechanisms = [
  *Mechanistic trade-offs, indication-specific sensitivity, and diminishing returns governing reference-guided bulk ICB deconvolution.*
  *(a)* Collinearity-accuracy phase portrait plotting mean bulk cohort response ROC-AUC against the log10 signature matrix condition number ($log_{10} kappa$) across all 312 reference configurations colored by Leiden clustering resolution, delineating the safe operating corridor ($kappa <= 40$, $"AUC" >= 0.70$) before severe multi-collinearity degrades regression stability.
  *(b)* Indication heterogeneity sensitivity matrix reporting Pearson correlations ($r$) between the 4 primary pipeline hyperparameters (cell states $K$, collinearity $kappa$, dataset count $K$, and Leiden resolution) and out-of-sample response ROC-AUC across each of the 9 bulk clinical trial validation cohorts.
  *(c)* Classifier model inductive bias $times$ reference granularity interaction across L2-penalized logistic regression, Random Forest, and Ridge linear probability models (LPM) from $K = 29$ to $K = 355$ cell states, demonstrating that $K = 42$ maximizes accuracy across all classifier architectures.
  *(d)* Single-cell subsampling titration curve ($18,000$ to $90,000$ cells) demonstrating logarithmic performance saturation beyond $approx 54,000$ profiled cells.
  *(e)* First-derivative marginal gain trajectory ($Delta "AUC" / Delta "res"$) across consecutive Leiden resolution increments ($0.25 arrow.r 2.00$), identifying the onset of diminishing returns above $"res" = 1.25$ and net accuracy loss at $"res" = 2.00$ ($-0.002$).
]

#let caption_supp_meta42_construction = [
  *Computational architecture and hierarchical Ward linkage assembly of the `Consensus-Meta42` pan-cancer reference atlas.*
  *(a)* Five-stage mathematical assembly flowchart illustrating the progression from 9 independent single-cell ICB cohorts ($N = 355$ local clusters across $G = 4,529$ union genes) through within-cohort relative specificity normalization, cross-cohort cosine similarity kernel estimation ($bold(S) in RR^(355 times 355)$), minimum-variance Ward hierarchical agglomerative fusion at $K^* = 42$, and final export of the well-conditioned reference matrix $bold(Phi)_("meta") in RR^(42 times 4,529)$ ($kappa = 18.4$).
  *(b)* Empirical $42 times 42$ cosine similarity matrix of the `Consensus-Meta42` states ordered by canonical lineage and internal Ward linkage, confirming sharp block-diagonal transcriptomic modularity.
  *(c)* Multi-cohort cluster fusion distribution quantifying how the 355 local single-cell clusters from the 9 source scRNA-seq studies coalesce into the 6 canonical tumor microenvironment lineages.
  *(d)* Hierarchical dendrogram cut optimization curve across $K in [10, 100]$ meta-states simultaneously evaluating out-of-sample bulk ICB response classification ROC-AUC and cross-modality single-cell/bulk concordance ($r$), establishing $K^* = 42$ as the Pareto-optimal granularity threshold.
]

#let caption_supp_meta42_atlas = [
  *Comprehensive transcriptomic, polar, and cross-modality characterization of the `Consensus-Meta42` reference atlas.*
  *(a)* Radial `Meta42` wheel displaying all 42 consensus cell states ordered clockwise by the 6 canonical TME lineages; radial spokes extend outward from the dashed zero-effect baseline ring ($hat(beta)_("meta") = 0$) for responder-enriched states ($hat(beta)_("meta") > 0$) and inward for non-responder barrier states ($hat(beta)_("meta") < 0$), with terminal node area proportional to single-cell Milo neighborhood enrichment magnitude ($|"logFC"|$) and perimeter callouts indicating representative state-marker pairs.
  *(b)* Cartesian block-diagonal reference signature matrix ($bold(Phi)_("meta")$) displaying standardized log2 expression specificity ($Z$-score) of top canonical HGNC marker genes across all 42 meta-states.
  *(c)* Cross-state co-regulation and antagonism landscape mapping single-cell Milo neighborhood logFC against bulk ICB responder effect ($hat(beta)_("meta")$) with transcriptomic coupling edges ($r > 0.48$), contrasting the responder effector-TLS hub (`M13` _CXCL13_+ Tex, `M14` _CXCL9_+ M1-TAM, `M15` _TCF7_+ stem CD8+, `M37` _MS4A1_+ TLS B, `M34` _GZMK_+ Tem) against the immunosuppressive stromal-myeloid barrier hub (`M23` _FAP_+ CAF, `M21` _FOXP3_+ Treg, `M19` _S100A8_+ MDSC).
  *(d)* Pan-cancer cohort ubiquity matrix ($42 times 9$) quantifying the number of fused local single-cell clusters contributing to each meta-state across the 9 source scRNA-seq cohorts.
]
