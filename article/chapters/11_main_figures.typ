// Main Article Figures (Figures 1–5; 1 figure = 1 dedicated page, placed after References)

#import "../captions/combined_markers_radar_triad_main_caption.typ": caption_combined_markers_radar_triad_main
#import "../captions/iatlas_epistasis_and_stromal_gating_main_caption.typ": caption_epistasis_stromal_gating_main
#import "../captions/sc9_reference_benchmarking_captions.typ": caption_sc9_summary_main

#set figure(placement: none)

= Figures

#figure(
  image("../figures/dataset_overview/figure_study_workflow_and_cohorts.svg", width: 100%, height: 185mm, fit: "contain"),
  caption: [
    *Overview of the tumor microenvironment ecosystem, clinical specimen processing, deconvolution and classification pipeline, and pan-cancer bulk and single-cell ICB cohorts.*
    *(a)* Cellular compartments of the heterogeneous tumor microenvironment (`TME`): (1) proliferating and hypoxic malignant cells, (2) effector immune populations (`CD8+` T, NK, and B lymphocytes), (3) suppressive immune subsets (`Tregs`, `M2` tumor-associated macrophages, and `MDSCs`), and (4) stromal fibroblasts (`CAFs`), extracellular matrix (`ECM`), and vasculature.
    *(b)* Clinical specimen and high-throughput sequencing workflow from surgical resection, pathological grossing, and spatial sampling through snap-frozen tissue preservation, RNA extraction, library preparation, and high-throughput sequencing.
    *(c)* Dual-stream single-cell and bulk deconvolution pipeline: single-cell transcriptomic atlases undergo quality control (`QC`), highly variable gene (`HVG`) selection, and principal-component harmonization (`PCA`/`Harmony`) to construct cell-state reference signature matrices ($Phi$), which are integrated with patient-cohort bulk RNA-seq profiles via Bayesian deconvolution (`BayesPrism`/`InstaPrism`) to identify driver cell states and train ensemble response classifiers.
    *(d)* Summary of the nine bulk RNA-seq ICB clinical trial cohorts ($n = 1,097$ evaluable patients across five solid cancer types; $357$ `RECIST v1.1` responders [`R`, $32.5%$] and $740$ non-responders [`NR`, $67.5%$]).
    *(e)* Summary of the nine single-cell ICB response cohorts in the `Milo` differential abundance compendium ($n = 613$ patient evaluations, $3,597,092$ cells across six solid tumor types; $304$ responders [$49.6%$] and $309$ non-responders [$50.4%$]).
  ]
) <fig-study-overview>

#pagebreak()

#figure(
  image("../figures/sc9-milo-deconv/fig_sc9_resolution_states_datasets_summary.svg", width: 100%),
  caption: caption_sc9_summary_main
) <fig-sc9-reference-summary>

#pagebreak()

#figure(
  image("../figures/signatures/fig_pr_auc_heatmap_transcriptomic_pre_standard.svg", width: 100%, height: 185mm, fit: "contain"),
  caption: [
    *Complete transcriptomic PR-AUC matrix for pre-treatment standard ICB response.* Precision–recall AUC of 28 RNA-seq gene-expression predictors (the deconvolution inference model and 27 literature biomarkers or single genes) evaluated across 9 individual clinical cohorts, grouped by cancer type (clear cell renal cell carcinoma, melanoma, pancreatic ductal, triple-negative breast, urothelial bladder), and 3 combined standardized cohort pools (melanoma $n = 338$, renal cell carcinoma $n = 263$, pan-cancer $n = 1,015$). Colour encodes PR-AUC (0.1–0.8); asterisks mark nominal significance as annotated in the panel.
  ]
) <fig-prauc-matrix>

#pagebreak()

#figure(
  image("../figures/milo-deconv/combined_main_figure_principled_gmrf.svg", width: 100%),
  caption: [
    *Cross-modality concordance of single-cell `milopy` differential abundance and bulk `InstaPrism` deconvolution, focused `Consensus-Meta42` poles, cross-cohort single-cell differential abundance, and Harmony-integrated single-cell embedding.*
    *(a)* `UMAP` of the outcome-blind `harmonypy` integration `I2` (cohort + sub-study) across all nine clinical ICB single-cell cohorts, coloured by cohort, with compartment centroids labelled (cohort `iLISI` $= 2.50$, compartment `cLISI` $= 1.00$, median of 3 cohorts per neighborhood); the full integration benchmark is shown in Supplementary Figure S35.
    *(b)* `UMAP` of the same `Harmony I2` single-cell atlas coloured by continuous cell-level `milopy` $log_2$ fold change (responder vs. non-responder; quasi-likelihood negative-binomial GLM `~ dataset_id + clinical_response` projected onto anchor cells), with compartment centroids labelled.
    *(c)* Single scatter of within-cohort $z$-score normalized single-cell `milopy` $overline("logFC")$ against within-cohort $z$-score normalized bulk CLR responder $hat(beta)$ (`Pan-iAtlas-Meta`) across all 9 cohort-specific single-cell ICB references (`Sade-Feldman`, `Gondal`, `Marteau`, `Yan`, `Hu`, `Liu`, `Hong`, `Luoma`, and `Raufi`; 42 `Consensus-Meta42` states each, colour = cohort) with a single combined OLS linear regression fit across all datasets; the per-cohort $3 times 3$ grid is shown in Supplementary Figure S57.
    *(d)* Paired cross-modality linear (`Pearson r`) and rank (`Spearman` $rho$) concordance summary across all 13 hierarchical references (including `Mel+PanCa`, `NSCLC-Uni`, `Unified-All9`, and `Meta-42`; all $p <= 5 times 10^(-4)$ across $B = 2,000$ permutations).
    *(e)* Merged cross-cohort `Milo` differential abundance summary across the 9 single-cell cohorts sharing the cohort axis: share of kNN neighborhoods significantly enriched in responders (positive) or in non-responders (negative) at $"FDR" < 0.10$ (left; 32,680 neighborhoods in total, of which 5,261 are responder-enriched and 2,439 non-responder-enriched) alongside the percentage of analysed cells lying in responder-enriched, non-responder-enriched, or non-significant `Milo` neighborhoods (right).
    *(f)* Coupled `Consensus-Meta42` extreme poles summary: dual-scale random-effects forest plot of the top 7 responder-enriched and top 7 non-responder-enriched `Consensus-Meta42` states ranked by `Pan-iAtlas-Meta` bulk deconvolution $hat(beta)_("meta")$ (95% CI bars, bottom axis) alongside single-cell `milopy` $overline("logFC")$ (top axis, 100% directional agreement across all 14 states; left) paired with the within-trial responder vs. non-responder CLR compositional shift heatmap ($hat(beta)_("deconv", b, k)$) across all 11 individual clinical bulk cohorts (`Hugo`, `Riaz`, `Liu`, `Gide`, `Rosenberg`, `Padron`, `Anders`, `McDermott`, `Choueiri`, `Jung`, and `Cho`; right).
  ]
) <fig-main-combined>

#pagebreak()

#figure(
  image("../figures/milo-deconv/combined_markers_radar_coulton.svg", width: 100%, height: 170mm, fit: "contain"),
  caption: caption_combined_markers_radar_triad_main
) <fig-metastate-deg-14>

#pagebreak()

#figure(
  image("../figures/milo-deconv/iatlas_epistasis_and_stromal_gating_main.svg", width: 100%, height: 170mm, fit: "contain"),
  caption: caption_epistasis_stromal_gating_main
) <fig-epistasis-stromal-gating-main>
