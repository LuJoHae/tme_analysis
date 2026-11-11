// Supplementary Figures (1 figure = 1 dedicated page, placed after Main Article Figures)

#import "../captions/iatlas_epistasis_supplementary_captions.typ": (
  caption_supp_pure_clr_all42,
  caption_supp_epistasis_power_benchmark,
  caption_supp_zero_leakage_gating_comparison,
  caption_supp_oof_diagnostics,
)
#import "../captions/coulton_tam_supplementary_captions.typ": (
  caption_supp_coulton_myeloid_similarity_ici,
  caption_supp_coulton_nearest_forest,
  caption_supp_coulton_all42_pan_lineage,
)
#import "../captions/tcell_caf_atlas_supplementary_captions.typ": (
  caption_supp_atlas_zheng_tcell,
  caption_supp_atlas_chu_tcell,
  caption_supp_atlas_shi_tnk,
  caption_supp_atlas_gao_caf,
  caption_supp_atlas_cords_caf,
  caption_supp_atlas_buechler_caf,
  caption_supp_atlas_shi_caf,
  caption_supp_triad_all42_specificity,
)
#import "../captions/meta42_radar_main_captions.typ": (
  caption_meta42_responder_radar,
  caption_meta42_patient_exemplars_radar,
)
#import "../captions/pan_cohort_pairwise_captions.typ": (
  caption_pan_cohort_core_hubs_main,
  caption_pan_cohort_unified_umap_main,
  caption_pan_cohort_pair_grid_main,
  caption_supp_pan_cohort_harmony_audit,
  caption_supp_pan_cohort_full_heatmap,
  caption_supp_pan_cohort_pair_umap,
)
#import "../captions/sc9_reference_benchmarking_captions.typ": (
  caption_supp_sc9_fig1,
  caption_supp_sc9_fig2,
  caption_supp_sc9_fig3,
  caption_supp_sc9_fig4,
  caption_supp_sc9_fig5,
)
#import "../captions/meta42_pipeline_supplementary_captions.typ": (
  caption_supp_pipeline_envelope,
  caption_supp_parameter_mechanisms,
  caption_supp_meta42_construction,
  caption_supp_meta42_atlas,
)
#import "../captions/deconv_approaches_comparison_supplementary_captions.typ": (
  caption_supp_deconv_factorial_matrix,
  caption_supp_deconv_validation_suite,
  caption_supp_deconv_head_to_head,
  caption_supp_deconv_gene_mechanisms,
)

#set figure(
  supplement: [Supplementary Figure],
  numbering: n => "S" + str(n),
  placement: none,
)
#counter(figure.where(kind: image)).update(0)

= Supplementary Figures

// S1–S3: Study Compendium, Milo Architecture & tme_datasets Pipeline DAG
#figure(
  image("../figures/dataset_overview/figure1_study_cohorts_overview.png", width: 100%, height: 162mm, fit: "contain"),
  caption: [
    *Comprehensive multi-omic immunotherapy dataset compendium, single-cell references, and analytical deconvolution framework.*
    *(a)* Multi-cohort dataset ecosystem comprising a 40,002-cell pan-cancer scRNA-seq reference (6 major lineages, 22 cell types, 58 fine-grained sub-clusters) alongside the matched clinical validation cohort (Sade-Feldman et al., $n=51$ pre/post ICB melanoma biopsies, 16,291 $"CD45"^+$ cells), 9 clinical trial cohorts ($n = 1,097$ patients) across melanoma, bladder, renal, pancreatic, and breast cancers, and 5 untreated primary tumor baseline cohorts from TCGA ($n = 2,932$ bulk tumors) defining $K=100$ unperturbed microenvironmental compositional centroids.
    *(b)* Multi-omic profiling layers spanning single-cell and bulk transcriptomics, whole-exome somatic mutations, subclonal variant allele frequency (VAF) spectra, tumor mutational burden (TMB), predictive response drivers (such as _TGM6_ alterations), standardized RECIST v1.1 clinical response criteria, and longitudinal treatment timing.
    *(c)* Unified analytical architecture executing single-cell continuous neighborhood testing (Milo kNN graph GLM), high-resolution pseudobulk deconvolution (`instaprism`), an automated multi-replicate distortion simulation suite ($N=5$, 485 runs across 8 single-variable modes, 4 compound regimes, and 2D interaction surfaces), and subclonal VAF joint optimization.
    *(d)* Key biomarker and methodological synthesis, illustrating empirical clinical concordance in matched melanoma ($rho = 0.853$), identification of concordant responder (naive B cells) and non-responder (macrophages) biomarkers alongside collinear cytotoxic T-cell discordance, primary deconvolution failure modes (patient marker dysregulation and cell size asymmetry), and quantification of the metastatic composition generalization gap.
  ]
) <fig-supp-study-overview>

#pagebreak()

#figure(
  image("../figures/deconvolution/how_milopy_works.png", width: 100%, height: 185mm, fit: "contain"),
  caption: [Methodological architecture of Milo / milopy differential abundance testing on single-cell k-nearest neighbor graphs. (1) kNN graph construction in PCA latent space and index vertex sampling with medoid refinement, bypassing discrete clustering boundaries. (2) Overlapping neighborhood definition, sparse binary incidence matrix $bold(M) in {0, 1}^(C times V)$, pairwise graph distance, and cross-cohort neighborhood composition $bold(C) = bold(B)^T bold(A)$. (3) Aggregation of single cells into patient-level count table $bold(N) in bb(N)_0^(V times J)$ with TMM library size normalization, eliminating single-cell pseudo-replication through hierarchical patient replicate modeling. (4) Quasi-likelihood Negative Binomial GLM with Empirical Bayes dispersion shrinkage, Spatial FDR p-value weighting via local graph connectivity density, and dual outputs (volcano plot and single-cell UMAP $log_2"FC"$ projection).]
) <fig-how-milopy-works>

#pagebreak()

#figure(
  image("../figures/dataset_overview/figure_tme_datasets_pipeline_dag.png", width: 100%, height: 158mm, fit: "contain"),
  caption: [
    *Unified `tme_datasets` registry queries and multi-stage computational provenance directed acyclic graph (DAG).*
    *(a) Tier I — Declarative Registry, Providers & Gene Panels (`packages/tme_datasets`):* Three foundational query domains expose (i) 9 single-cell ICB response cohorts (`SC-9`, $3.60 times 10^6$ cells across melanoma, pan-cancer, colorectal, triple-negative breast, non-small cell lung, head-and-neck squamous, and pancreatic ductal cancers) via `load_dataset()`, `resolve_milopy_cohorts()`, and `load_backed()`; (ii) 11 bulk RNA-seq clinical trials plus 4 indication/pan-cancer pooled strata (`BULK-15`, $N = 1,058$ evaluable biopsies across 9 iAtlas trials and 2 NSCLC cohorts) via `load_dataset()` and `load_iatlas_cohort_or_combined()`; and (iii) 17+ auxiliary single-cell reference atlases, 5 untreated TCGA primary tumor baselines ($N = 2,932$), and curated gene collections (`IMMUNOTHERAPY_GENE_PANEL`, `IMMUNE_CHECKPOINT_GENES`, Bagaev TME signatures) via `load_geneset_collection()` and `list_registered_datasets()`.
    *(b) Tier II — Primary Query, QC & `.h5ad` Serialization Scripts (`scripts/manuscript_pipeline/`):* `1a_download_and_preprocess_all.py`, `1b_batch_preprocess_cohorts.py`, and `1c_run_qc_pipeline_all_single_cell.py` `[1A–1C]` execute Ensembl v111 symbol reconciliation and write atomic sparse `AnnData` caches (`data/preprocessed/{dataset_id}.h5ad`) alongside per-script summaries under `output/manuscript_pipeline/<script_stem>/`; in parallel, `2a1_sample_and_analyze_cohorts.py` and `2a2_combine_cohorts_and_analyze.py` `[2A]` query `SC-9` to compute single-cell kNN Negative Binomial GLM differential abundance (`output/manuscript_pipeline/2a1_sample_and_analyze_cohorts/*.parquet`), while `2b_run_stability_selection_iatlas.py` `[2B]` queries `BULK-15` and the 101-gene IO panel across 12 cohorts $times$ 13 fitters (`output/manuscript_pipeline/2b_run_stability_selection_iatlas/*.parquet`).
    *(c) Tier III — Chained Analytical Compute Engines:* `3a_sade_feldman_pseudobulk_deconv_benchmark.py` `[3A]` constructs 355 Leiden ($r = 2.0$) reference centroids; `3b_pan_cohort_unified_deconv_benchmark.py` `[3B]` Ward-clusters those 355 states into the 42-state `Consensus-Meta42` pan-cancer atlas (`M01`–`M42`); `3c_iatlas_instaprism_milopy_benchmark.py` `[3C]` emits 16 Parquet tables across bulk linearization, `InstaPrism` EM deconvolution ($13 times 15$), 5-family ML response prediction, zero-shot `LOCO` transfer, `Pan-iAtlas-Meta` concordance, 168-gene marker specificity, calibrated unblended `Pure CLR` / `SC-Markers` epistasis, and zero-leakage `5`-fold `OOF` stromal gating; joined by `[3D]`–`[3G]`.
    *(d) Tier IV — Publication Figure, Caption & Manuscript Synthesis:* `4b_plot_iatlas_instaprism_milopy_figures.py` `[4B]` consumes all 16 Parquet tables from `[3C]` to render 14 heading-free and atlas figures plus standalone `.typ` caption modules (`article/captions/*.typ`), joined by `[4A]`, `[4C]`, `[4D]`, and `[4E]` (`4e_generate_tme_datasets_pipeline_dag.py`) in `article/Makefile` for 300 DPI `rsvg-convert` rasterization and Typst compilation into `article/main.pdf`.
  ]
) <fig-tme-datasets-pipeline-dag>

#pagebreak()

// S33–S39: Collinearity-Regularized Deconvolution & 3-Factorial Synthetic Benchmarks
#figure(
  image("../figures/deconvolution/synthetic_dataset_overview.png", width: 100%, height: 184mm, fit: "contain"),
  caption: [Controlled synthetic benchmark dataset architecture. (A) Reference signature correlation matrices across sibling correlation levels ($r = 0.0, 0.6, 0.9, 0.99$) across $S=12$ synthetic states and $T=4$ lineages (`Lineage_1` to `Lineage_4`), demonstrating emerging $3 times 3$ intra-lineage block collinearity while preserving cross-lineage orthogonality ($r approx 0.0$). (B) Ground truth non-uniform Dirichlet proportion distribution across the 12 synthetic subsets (`L1_State_1` through `L4_State_3`, $N=50$ samples, $bold(alpha) = [1.2, 1.0, 0.8, 1.0, 0.8, 0.6, 1.2, 1.0, 0.8, 0.8, 0.6, 0.4]$). (C) Representative bulk mixture sequencing read counts (Sample 01, depth $= 80,000$ reads across synthetic genes).]
) <fig-synthetic-datasets>

#pagebreak()

#figure(
  image("../figures/deconvolution/collinearity_regularized_synthetic_benchmark.png", width: 100%, height: 168mm, fit: "contain"),
  caption: [Continuous-spectrum statistical benchmark of collinearity-aware regularized deconvolution across sibling correlations ($r in [0.0, 0.99]$, $S=12$ states, $N=50$ samples per point) across all seven methods: unregularized NNLS, RegDeconv, Rectangle (DWLS-QP), CIBERSORT (reimplementation, $nu$-SVR), InstaPrism, BayesPrism (Gibbs), and the official CIBERSORTx container (`cibersortx/fractions:latest`). (A) Reference Hessian condition number $kappa(bold(H))$ on $log_(10)$ scale: standard Poisson likelihood escalates toward $3.8 times 10^3$ as $r -> 0.99$, whereas Adaptive Deficit Regularization bounds $kappa <= 2,000$ ($kappa = 1,618$). (B) State-level proportion Mean Squared Error (MSE, Mean $plus.minus$ SD across $N=50$ samples): at $r <= 0.8$, all methods maintain minimal error ($2.5 times 10^(-5)$ to $1.8 times 10^(-4)$); at $r = 0.99$, regularized and Bayesian methods contract collinear siblings toward their centroid, trading small shrinkage bias for topological stability. Dashed slate line confirms broad lineage MSE ($approx 5.4 times 10^(-4)$) is invariant across methods. (C) Intra-lineage sibling cross-talk correlation ($"Corr"(hat(theta)_a, hat(theta)_b)$): BayesPrism and InstaPrism scale smoothly toward $+1.00$, RegDeconv eliminates anti-correlation and jumps to $+0.44$ at $r = 0.99$, whereas NNLS ($-0.130$), Rectangle ($0.000$), CIBERSORT (reimpl., $-0.135$), and CIBERSORTx Docker ($-0.132$) remain trapped in negative cross-talk or zero correlation. (D) Spurious state dropout rate (% true non-zero states collapsed to boundary $0.0$): NNLS collapses $10.0%$, CIBERSORT reimpl. $9.5%$, and CIBERSORTx Docker $7.5%$ of states to false zeros, whereas RegDeconv suppresses dropouts to $0.5%$, and BayesPrism and InstaPrism maintain $0.0%$ dropouts.]
) <fig-collinearity-benchmark>

#pagebreak()

#figure(
  image("../figures/deconvolution/deconvolution_results_vs_ground_truth.png", width: 100%, height: 184mm, fit: "contain"),
  caption: [Deconvolution ground-truth recovery under severe collinearity ($r = 0.99$). (A) Scatter plot grid of inferred cell fraction ($hat(theta)$) vs. ground truth ($theta^*$) across all six methods with diagonal $y=x$ dashed reference lines, showing $R^2$, RMSE, and Pearson $r$. (B) True vs. predicted cellular composition for representative synthetic samples (Samples 01 to 06), illustrating how unregularized NNLS collapses sibling state mass to zero while RegDeconv, BayesPrism, and InstaPrism preserve proportional representation across synthetic subsets `L1_State_1` to `L4_State_3`.]
) <fig-ground-truth-comparison>

#pagebreak()

#figure(
  image("../figures/deconvolution/variance_and_dropout_diagnostics.png", width: 100%, height: 178mm, fit: "contain"),
  caption: [Estimator bias-variance decomposition and biological zero dropout stress-testing. (A1) Exact MSE decomposition ($"MSE" = "Bias"^2 + "Variance"$) across $B=30$ independent sequencing draws from the same biological mixture at $r = 0.99$: unregularized NNLS is $95.6%$ estimator variance, whereas BayesPrism and InstaPrism are $>99.5%$ centroid shrinkage bias with near-zero variance. (A2) Estimator variance $"Var"(hat(theta))$ on $log_(10)$ scale: BayesPrism ($2.8 times 10^(-5)$) and InstaPrism ($3.0 times 10^(-6)$) reduce technical noise variance by $45times$ to $400times$ relative to NNLS ($1.2 times 10^(-3)$). (B1) False positive phantom detection mass (%) under true biological zero absence ($theta^* = 0$): Dirichlet shrinkage priors leak $8.3%$ mass into absent sibling states when their lineage is active, and BayesPrism accumulates $15.7%$ total phantom mass across absent lineages. (B2) Reconstruction RMSE on active non-zero subsets.]
) <fig-variance-dropout>

#pagebreak()

#figure(
  image("../figures/deconvolution/bayesprism_superior_regimes.png", width: 100%, height: 170mm, fit: "contain"),
  caption: [Empirical regimes of Bayesian deconvolution dominance across biological and technical stress tests under severe collinearity ($r = 0.99$). (A) Low sequencing depth ($n_("total") in [1,000, 80,000]$): as read depth drops below $10"k"$, Poisson sampling noise detonates the variance of NNLS ($"MSE" = 0.0152$), Rectangle ($"MSE" = 0.0139$), and CIBERSORT ($"MSE" = 0.0133$), while InstaPrism ($"MSE" = 0.0045$) and BayesPrism ($"MSE" = 0.0066$) achieve up to $3.4times$ lower error. (B) Patient-specific reference expression drift ($sigma_("drift") in [0.0, 0.75]$): multiplicative lognormal distortion of tumor/stromal signatures causes linear least-squares error to surge to $0.0113$, whereas multinomial likelihood dampens gene expression outliers, keeping InstaPrism ($0.0048$) and BayesPrism ($0.0051$) completely flat ($2.4times$ error reduction). (C) Biological overdispersion and transcriptional bursting ($alpha_("disp") in [0.0, 0.30]$): under Negative Binomial extra-Poisson variance, $L_2$-loss penalties over-leverage bursty genes ($"MSE" = 0.0153$), whereas Bayesian count models bound outlier leverage ($"MSE" = 0.0052$, $3.0times$ reduction). (D) Balanced intra-lineage plasticity ($alpha_("within") in [1, 40]$): when collinear sibling states co-occur evenly around parity ($theta_(t,s) approx 1/3$), shrinkage bias drops to zero, and InstaPrism ($"MSE" = 0.00017$, $13.9times$ reduction vs. NNLS) and BayesPrism ($"MSE" = 0.00036$, $6.7times$ reduction) decisively dominate all linear models.]
) <fig-superior-regimes>

#pagebreak()

#figure(
  image("../figures/deconvolution/unified_bias_variance_landscape.png", width: 100%, height: 182mm, fit: "contain"),
  caption: [Two-dimensional faceted bias-variance decomposition landscape ($"MSE" = "Bias"^2 + "Variance"$) across collinearity levels $r in [0.0, 0.99]$ (columns) and sequencing count parameters $n_("total") in [2,500, 80,000]$ (rows) across all seven deconvolution methods under heterogeneous biological proportions. Stacked bars represent total state MSE, partitioned into estimator variance $"Var"(hat(theta))$ (purple) and squared bias $"Bias"^2(hat(theta))$ (gold). At low counts ($n_("total") = 2,500$) and severe collinearity ($r = 0.99$), unregularized NNLS, CIBERSORT, and CIBERSORTx Docker explode into pure estimator variance, whereas InstaPrism and BayesPrism maintain near-zero variance, achieving lower total MSE.]
) <fig-unified-landscape>

#pagebreak()

#figure(
  image("../figures/deconvolution/algorithm_dominance_matrix.png", width: 100%, height: 184mm, fit: "contain"),
  caption: [Empirical algorithm dominance decision map identifying the winning deconvolution method (lowest total MSE) across collinearity $r$ and count parameter $n_("total")$ for each biological state architecture. (Left) Balanced intra-lineage plasticity: InstaPrism dominates almost universally, with RegDeconv optimal at high counts and moderate collinearity. (Center) Heterogeneous proportions: InstaPrism dominates at low counts ($n_("total") <= 10,000$), while NNLS and RegDeconv share the high-count regime. (Right) Structural state dropout ($theta^* = 0$): RegDeconv dominates across moderate-to-high counts by avoiding Dirichlet phantom leakage, while InstaPrism buffers against Poisson shot noise at low counts.]
) <fig-dominance-map>

#pagebreak()

// S40–S49: Cross-Modality ICB Prediction, Epistasis Calibration, Stromal Gating & Coulton TAM Validation
#figure(
  image("../figures/milo-deconv/iatlas_response_prediction_auc_heatmap_and_bars.png", width: 100%, height: 182mm, fit: "contain"),
  caption: [
    *Multi-tier cross-modality ICB clinical response prediction (`ROC-AUC`) across 13 single-cell references and 15 bulk cohorts/strata ($N = 1,058$ patients).*
    *(a)* Cross-validated `ROC-AUC` heatmap ($13$ single-cell references spanning Tiers 1–3 $times$ $15$ individual clinical trials and pooled strata) evaluated with $ell_2$-penalized logistic regression on compositional CLR and $Delta$CLR log-ratios after exact bulk linearization ($2^X - 1 -> 10^6$ TPM).
    *(b)* Out-of-fold `ROC-AUC` comparison on the 42-state pan-cancer `Consensus-Meta42` reference across all 15 bulk cohorts for the three primary classifier paradigms (`Logistic Reg.`, `Linear Reg.`, and `Random Forest`; full 5-family regularized metrics including `Elastic-Net` and `Lasso L1` reported in text), peaking at $0.945$ in `Cho`, $0.782$ in `NSCLC Pool`, $0.769$ in `Hugo`, and $0.764$ in `Gide`.
  ]
) <fig-iatlas-roc-auc>

#pagebreak()

#figure(
  image("../figures/milo-deconv/iatlas_response_prediction_prauc_heatmap_and_bars.png", width: 100%, height: 182mm, fit: "contain"),
  caption: [
    *Prevalence-calibrated precision–recall performance (`PR-AUC` and excess gain $Delta "PR-AUC"$) across 13 references and 15 bulk cohorts.*
    *(a)* Cross-validated `PR-AUC` heatmap across all $13 times 15$ reference-cohort pairs (`Logistic L2`).
    *(b)* Multi-model out-of-fold `PR-AUC` comparison (`Logistic Reg.`, `Linear Reg.`, and `Random Forest`) across all 15 cohorts relative to each cohort's positive responder prevalence baseline $pi_("pos")$ (vertical tick markers), demonstrating consistent precision lift across class-imbalanced clinical trials (`PanCancer Pool` `PR-AUC` $= 0.423$ vs. baseline $pi_("pos") = 0.314$, $Delta "PR-AUC" = +0.109$; mean $Delta "PR-AUC" = +0.260$ across cohorts for `Linear Reg.` on `Consensus-Meta42` vs. $+0.191$ on `Unified-All9`).
  ]
) <fig-iatlas-prauc>

#pagebreak()

#figure(
  image("../figures/milo-deconv/iatlas_nonlinear_conditional_interactions_all42_pure_clr.png", width: 100%, height: 175mm, fit: "contain"),
  caption: caption_supp_pure_clr_all42
) <fig-supp-s19>

#pagebreak()

#figure(
  image("../figures/milo-deconv/iatlas_epistasis_discovery_power_benchmark.png", width: 100%, height: 172mm, fit: "contain"),
  caption: caption_supp_epistasis_power_benchmark
) <fig-supp-s20>

#pagebreak()

#figure(
  image("../figures/milo-deconv/iatlas_zero_leakage_gating_comparison.png", width: 100%, height: 172mm, fit: "contain"),
  caption: caption_supp_zero_leakage_gating_comparison
) <fig-supp-s21>

#pagebreak()

#figure(
  image("../figures/milo-deconv/iatlas_oof_supplementary_diagnostics.png", width: 100%, height: 175mm, fit: "contain"),
  caption: caption_supp_oof_diagnostics
) <fig-supp-s22>

#pagebreak()

#figure(
  image("../figures/coulton-tam/fig1_meta42_myeloid_coulton_similarity_and_ici.png", width: 100%, height: 172mm, fit: "contain"),
  caption: caption_supp_coulton_myeloid_similarity_ici
) <fig-supp-s23>

#pagebreak()

#figure(
  image("../figures/coulton-tam/fig2_macrophage_nearest_substate_ici_forest.png", width: 100%, height: 174mm, fit: "contain"),
  caption: caption_supp_coulton_nearest_forest
) <fig-supp-s24>

#pagebreak()

#figure(
  image("../figures/coulton-tam/fig3_all42_pan_lineage_coulton_specificity.png", width: 100%, height: 174mm, fit: "contain"),
  caption: caption_supp_coulton_all42_pan_lineage
) <fig-supp-s25>

#pagebreak()

// S50–S51: Finite-Sample Controlled Stability Selection & Multi-Paradigm Fitter Benchmark
#figure(
  image("../figures/stability-selection/all_cohorts_stability_paths_composite.png", width: 100%, height: 184mm, fit: "contain"),
  caption: [Finite-sample controlled stability paths across 12 clinical immunotherapy cohorts ($n=1,015$ total patients). Each panel depicts empirical feature selection probabilities $hat(Pi)_k^lambda$ across the budgeted regularization path for an individual cohort or multi-cohort pooled setting. The horizontal dashed red line marks the stability selection threshold ($pi_("thr") = 0.60$ for single cohorts, $pi_("thr") = 0.70$ for multi-cohort pooled regimes). Solid colored paths denote stable biomarkers surpassing $pi_("thr")$ under strict PFER error control ($"PFER" <= 2.5$), while faint gray paths denote unselected noise features. The universal chemokine marker _CXCL13_ (highlighted in gold/orange) consistently achieves maximum stability across indications.]
) <fig-stability-paths-composite>

#pagebreak()

#figure(
  image("../figures/stability-selection/fitter_gene_overlap_nature.png", width: 100%, height: 182mm, fit: "contain"),
  caption: [Cross-paradigm benchmark of 13 machine learning and econometric fitters in stability selection. (A) Gene selection frequency across all 13 fitters. _CXCL13_ emerged as the invariant, model-agnostic driver of response, selected by 12 out of 13 algorithms. Core immune effectors (_CD8A_, _CXCL9_, _CXCL10_, _IFNG_, _LAG3_, _TIGIT_) were recurrently recovered across multiple paradigms. (B) Fitter paradigm classification highlighting convex GLMs (Lasso, ElasticNet), non-convex/adaptive estimators (SCAD, MCP, AdaptiveLasso, HardThresholding), tree ensembles (LightGBM, XGBoost, RandomForest), cohort-aware econometric models (LassoFWL, MERF), and ordered weighted $ell_1$ estimators (OSCAR, SLOPE). (C) Gene selection matrix across all fitters, demonstrating consensus patterns and method-specific inductive biases.]
) <fig-fitter-overlap>

#pagebreak()

// S52–S56: Multi-Cohort Reference Scaling, Granularity Regimes & Classifier Selection Framework
#figure(
  image("../figures/sc9-milo-deconv/fig1_predictive_capacity_benchmark.png", width: 100%, height: 184mm, fit: "contain"),
  caption: caption_supp_sc9_fig1
) <fig-supp-sc9-fig1>

#pagebreak()

#figure(
  image("../figures/sc9-milo-deconv/fig2_signature_collinearity_benchmark.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_sc9_fig2
) <fig-supp-sc9-fig2>

#pagebreak()

#figure(
  image("../figures/sc9-milo-deconv/fig3_reference_selection_decision_framework.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_sc9_fig3
) <fig-supp-sc9-fig3>

#pagebreak()

#figure(
  image("../figures/sc9-milo-deconv/fig4_winning_classifier_biomarker_drivers.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_sc9_fig4
) <fig-supp-sc9-fig4>

#pagebreak()

#figure(
  image("../figures/sc9-milo-deconv/fig5_selected_cohort_auc_and_marginal_gain.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_sc9_fig5
) <fig-supp-sc9-fig5>

#pagebreak()

#figure(
  image("../figures/meta42-pipeline/fig_a_parameter_optimization_envelope.png", width: 100%, height: 178mm, fit: "contain"),
  caption: caption_supp_pipeline_envelope
) <fig-supp-meta42-envelope>

#pagebreak()

#figure(
  image("../figures/meta42-pipeline/fig_b_parameter_mechanisms_and_tradeoffs.png", width: 100%, height: 178mm, fit: "contain"),
  caption: caption_supp_parameter_mechanisms
) <fig-supp-meta42-mechanisms>

#pagebreak()

#figure(
  image("../figures/meta42-pipeline/fig_c_meta42_computational_construction.png", width: 100%, height: 178mm, fit: "contain"),
  caption: caption_supp_meta42_construction
) <fig-supp-meta42-construction>

#pagebreak()

#figure(
  image("../figures/meta42-pipeline/fig_d_meta42_atlas_polar_and_signatures.png", width: 100%, height: 178mm, fit: "contain"),
  caption: caption_supp_meta42_atlas
) <fig-supp-meta42-atlas>

// External Pan-Cancer T-Cell, CAF, and Unified 42-State Triad Reference Atlas Validation

#pagebreak()

#figure(
  image("../figures/pan-cancer-atlases/meta42_Zheng_TCell_Atlas_similarity_and_ici.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_atlas_zheng_tcell
) <fig-supp-atlas-zheng-tcell>

#pagebreak()

#figure(
  image("../figures/pan-cancer-atlases/meta42_Chu_TCell_Atlas_similarity_and_ici.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_atlas_chu_tcell
) <fig-supp-atlas-chu-tcell>

#pagebreak()

#figure(
  image("../figures/pan-cancer-atlases/meta42_Shi_PanCancer_Atlas_TNK_similarity_and_ici.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_atlas_shi_tnk
) <fig-supp-atlas-shi-tnk>

#pagebreak()

#figure(
  image("../figures/pan-cancer-atlases/meta42_Gao_CAF_Atlas_similarity_and_ici.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_atlas_gao_caf
) <fig-supp-atlas-gao-caf>

#pagebreak()

#figure(
  image("../figures/pan-cancer-atlases/meta42_Cords_CAF_Atlas_similarity_and_ici.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_atlas_cords_caf
) <fig-supp-atlas-cords-caf>

#pagebreak()

#figure(
  image("../figures/pan-cancer-atlases/meta42_Buechler_CAF_Atlas_similarity_and_ici.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_atlas_buechler_caf
) <fig-supp-atlas-buechler-caf>

#pagebreak()

#figure(
  image("../figures/pan-cancer-atlases/meta42_Shi_PanCancer_Atlas_CAF_similarity_and_ici.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_atlas_shi_caf
) <fig-supp-atlas-shi-caf>

#pagebreak()

#figure(
  image("../figures/pan-cancer-atlases/meta42_triad_all42_external_specificity.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_triad_all42_specificity
) <fig-supp-triad-all42>

// Pan-cohort Harmony integration audit, full cross-compartment heatmap & broad cell-type pair UMAPs

#pagebreak()

#figure(
  image("../figures/milo-pairwise/pan_cohort_harmony_integration_audit.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_pan_cohort_harmony_audit
) <fig-supp-pan-harmony-audit>

#pagebreak()

#figure(
  image("../figures/milo-pairwise/pan_cohort_cross_compartment_z_acat_heatmap.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_pan_cohort_full_heatmap
) <fig-supp-pan-full-heatmap>

#pagebreak()

#figure(
  image("../figures/milo-pairwise/pan_cohort_umap_t_nk_vs_caf.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_pan_cohort_pair_umap([`T-NK` $<->$ `CAF`], 6, 2)
) <fig-supp-pan-umap-t-nk-vs-caf>

#pagebreak()

#figure(
  image("../figures/milo-pairwise/pan_cohort_umap_t_nk_vs_endothelial.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_pan_cohort_pair_umap([`T-NK` $<->$ `Endothelial`], 2, 2)
) <fig-supp-pan-umap-t-nk-vs-endothelial>

#pagebreak()

#figure(
  image("../figures/milo-pairwise/pan_cohort_umap_myeloid_macrophages_vs_caf.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_pan_cohort_pair_umap([`Myeloid/Macrophages` $<->$ `CAF`], 2, 2)
) <fig-supp-pan-umap-myeloid-macrophages-vs-caf>

#pagebreak()

#figure(
  image("../figures/milo-pairwise/pan_cohort_umap_myeloid_macrophages_vs_endothelial.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_pan_cohort_pair_umap([`Myeloid/Macrophages` $<->$ `Endothelial`], 2, 2)
) <fig-supp-pan-umap-myeloid-macrophages-vs-endothelial>

#pagebreak()

#figure(
  image("../figures/milo-pairwise/pan_cohort_umap_b_cells_vs_caf.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_pan_cohort_pair_umap([`B cell` $<->$ `CAF`], 1, 0)
) <fig-supp-pan-umap-b-cells-vs-caf>

#pagebreak()

#figure(
  image("../figures/milo-pairwise/pan_cohort_umap_b_cells_vs_endothelial.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_pan_cohort_pair_umap([`B cell` $<->$ `Endothelial`], 0, 1)
) <fig-supp-pan-umap-b-cells-vs-endothelial>

#pagebreak()

#figure(
  image("../figures/milo-pairwise/pan_cohort_umap_t_nk_vs_myeloid_macrophages.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_pan_cohort_pair_umap([`T-NK` $<->$ `Myeloid/Macrophages`], 4, 14)
) <fig-supp-pan-umap-t-nk-vs-myeloid-macrophages>

#pagebreak()

#figure(
  image("../figures/milo-pairwise/pan_cohort_umap_b_cells_vs_myeloid_macrophages.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_pan_cohort_pair_umap([`B cell` $<->$ `Myeloid/Macrophages`], 2, 0)
) <fig-supp-pan-umap-b-cells-vs-myeloid-macrophages>

#pagebreak()

#figure(
  image("../figures/milo-pairwise/pan_cohort_umap_t_nk_vs_b_cells.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_pan_cohort_pair_umap([`T-NK` $<->$ `B cell`], 1, 2)
) <fig-supp-pan-umap-t-nk-vs-b-cells>

#pagebreak()

#figure(
  image("../figures/milo-pairwise/pan_cohort_umap_t_nk_vs_plasma.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_pan_cohort_pair_umap([`T-NK` $<->$ `Plasma cell`], 1, 1)
) <fig-supp-pan-umap-t-nk-vs-plasma>

#pagebreak()

#figure(
  image("../figures/milo-pairwise/pan_cohort_umap_b_cells_vs_plasma.png", width: 100%, height: 180mm, fit: "contain"),
  caption: caption_supp_pan_cohort_pair_umap([`B cell` $<->$ `Plasma cell`], 1, 0)
) <fig-supp-pan-umap-b-cells-vs-plasma>

#pagebreak()

#figure(
  image("../figures/milo-deconv/iatlas_cross_cohort_concordance_transfer_matrix.png", width: 100%, height: 182mm, fit: "contain"),
  caption: [
    *Cross-indication effect-size concordance matrix ($13 times 16$) and zero-shot Leave-One-Cohort-Out (`LOCO`) clinical transferability.*
    *(a)* Complete $13$ (single-cell reference) $times$ $16$ (`Pan-iAtlas Meta (9)`, 4 indication/pan-cancer pools, and 11 individual bulk trials) linear Pearson $r$ concordance matrix between single-cell `milopy` $overline("logFC")_k$ and bulk `InstaPrism` CLR $hat(beta)_("deconv", k)$ across all 208 reference-target pairs.
    *(b)* Within-cohort $5 times 5$ CV vs. zero-shot Leave-One-Cohort-Out (`LOCO`) cross-trial `ROC-AUC` generalization across all 11 individual clinical trials (trained strictly on the other 10 trials without cohort refitting; `Random Forest` mean `LOCO ROC-AUC` $= 0.661$, peaking at $0.945$ in `Cho`, $0.779$ in `Gide`, $0.744$ in `Choueiri`, and $0.730$ in `Jung`).
  ]
) <fig-iatlas-transfer-loco>

#pagebreak()

#figure(
  image("../figures/milo-deconv/iatlas_consensus_metastate_all42_forest_and_tme_profile.png", width: 100%, height: 178mm, fit: "contain"),
  caption: [
    *Complete 42-state `Consensus-Meta42` pan-cancer atlas (`M01`–`M42`): full-spectrum cross-modality forest plot and 11-trial TME profile.*
    *(a)* Exhaustive dual-scale random-effects forest plot across all 42 `Consensus-Meta42` meta-states ordered monotonically from top responder-enriched (`M13: CXCL13+PDCD1+ Tex`, $hat(beta)_("meta") = +2.70$) through intermediate plastic/homeostatic states (`M01`–`M42` continuum; 26 states with $hat(beta)_("meta") > 0$, 16 states with $hat(beta)_("meta") < 0$) down to top non-responder-enriched (`M23: FAP+TGFB1+ Mat. CAF`, $hat(beta)_("meta") = -3.13$), overlaying bulk `Pan-iAtlas-Meta` $hat(beta)_("meta")$ (95% CIs, bottom axis) with single-cell `milopy` $overline("logFC")$ (top axis; $r = +0.928, rho = +0.917$, $88.1%$ sign concordance across 37/42 states).
    *(b)* Complete $42 times 11$ within-trial responder vs. non-responder CLR shift matrix across all 11 individual clinical trials, revealing conserved cross-indication polarity alongside lineage-specific modulation in renal cell (`McDermott`, `Choueiri`) and pancreatic (`Padron`) cohorts.
  ]
) <fig-metastate-forest-all42>

#pagebreak()

#figure(
  image("../figures/milo-deconv/iatlas_consensus_metastate_all42_top_deg_markers.png", width: 100%, height: 182mm, fit: "contain"),
  caption: [
    *Complete 42-state `Consensus-Meta42` marker gene compendium (`168` protein-coding HGNC state-gene pairs across `M01`–`M42`).*
    *(a)* Exhaustive $42 times 4$ HGNC marker specificity matrix (`Rank #1`–`Rank #4` in $bold(Phi)_"meta"$) ordered monotonically by `Pan-iAtlas-Meta` $hat(beta)_("meta")$ across all 26 responder-leaning states ($hat(beta)_("meta") > 0$) and all 16 non-responder-leaning states ($hat(beta)_("meta") < 0$).
    *(b)* Pan-cancer bulk RNA-seq 4-gene signature responder vs. non-responder Welch $t$-statistics ($N = 1,015$ patients) for every meta-state, demonstrating direct alignment with independent gene-level stability selection (Sections 3.22–3.25), where _CXCL13_ (`M13`, `M24`), _CXCL9_/_STAT1_ (`M14`), and cytotoxic effectors emerge as invariant bulk biomarkers.
  ]
) <fig-metastate-deg-all42>

#pagebreak()

#figure(
  image("../figures/milo-deconv/iatlas_consensus_metastate_responder_radar.png", width: 100%, height: 178mm, fit: "contain"),
  caption: caption_meta42_responder_radar
) <fig-meta42-responder-radar>

#pagebreak()

#figure(
  image("../figures/milo-deconv/iatlas_consensus_metastate_patient_exemplars_radar.png", width: 100%, height: 178mm, fit: "contain"),
  caption: caption_meta42_patient_exemplars_radar
) <fig-meta42-patient-exemplars-radar>

#pagebreak()

#figure(
  image("../figures/signatures/fig_roc_auc_heatmap_transcriptomic_pre_standard.svg", width: 100%, height: 185mm, fit: "contain"),
  caption: [
    *Complete transcriptomic ROC-AUC matrix for pre-treatment standard ICB response.* Area under the ROC curve of the same 28 RNA-seq gene-expression predictors (the deconvolution inference model and 27 literature biomarkers or single genes) across the 9 individual clinical cohorts, grouped by cancer type, and the 3 combined standardized cohort pools (melanoma $n = 338$, renal cell carcinoma $n = 263$, pan-cancer $n = 1,015$). Colour encodes ROC-AUC (diverging scale centred near 0.6); asterisks mark nominal significance as annotated in the panel. This is the ROC-AUC counterpart of the PR-AUC matrix in @fig-prauc-matrix.
  ]
) <fig-supp-roc-auc-matrix>

#pagebreak()

#figure(
  image("../figures/dataset_overview/figure_datasets_tabular_compendium.svg", width: 100%, height: 185mm, fit: "contain"),
  caption: [
    *Single-cell immunotherapy response compendium: cohorts, indications, and clinical therapies.* Tabular summary of the 9 single-cell ICB cohorts (6 indications; 613 patients, 304 responders and 309 non-responders; 371,568 analyzed cells), listing for each cohort the publication, cancer type, treatment and target mechanism, sequencing platform and chemistry, biopsy timing, analyzed cell number, RECIST 1.1 response counts, and objective response rate (ORR).
  ]
) <fig-supp-sc-compendium-table>

#pagebreak()

#figure(
  image("../figures/single-cell/milopy_cross_cohort_volcano_grid.png", width: 100%, height: 185mm, fit: "contain"),
  caption: [
    *Single-cell differential abundance testing across 9 immunotherapy response cohorts.* Volcano plots (log#sub[2] fold change, responder vs. non-responder, against $-log_10("FDR")$) of `milopy` quasi-likelihood negative-binomial GLM tests on kNN neighborhoods for each cohort (a: GSE120575, melanoma; b: CELLxGENE_7b20c613, melanoma; c: CELLxGENE_05a8c945, colorectal; d: CELLxGENE_6f9de485, breast; e: GSE207422, NSCLC; f: GSE243013, NSCLC; g: GSE233203, NSCLC; h: GSE200996, HNSCC; i: GSE316195, PDAC). The horizontal dashed rule marks FDR = 0.10. Neighborhoods are coloured as enriched in responders, enriched in non-responders, or not significant; per-panel counts of up- and down-regulated neighborhoods are given in the panel headers.
  ]
) <fig-supp-milo-volcano-grid>

#pagebreak()

#figure(
  image("../figures/single-cell/milopy_cross_cohort_umaps_logfc_bitmap.svg", width: 100%, height: 185mm, fit: "contain"),
  caption: [
    *Single-cell landscape of immunotherapy response across 9 solid tumor cohorts.* Two-dimensional UMAP embeddings of each cohort (panels a–i as in the preceding volcano figure; 15,142–50,000 cells per cohort after subsampling), coloured by the continuous inferred `milopy` log#sub[2] fold change (responder vs. non-responder, clamped at $plus.minus 5$).
  ]
) <fig-supp-milo-umaps-logfc>

#pagebreak()

#figure(
  image("../figures/milo-deconv/iatlas_milopy_vs_deconv_scatter_grid_supp_principled_gmrf.svg", width: 100%, height: 185mm, fit: "contain"),
  caption: [
    *Per-cohort cross-modality concordance between single-cell `milopy` differential abundance and bulk `InstaPrism` deconvolution effect sizes.* $3 times 3$ scatter grid across all 9 cohort-specific single-cell ICB references (a: `Sade-Feldman`, b: `Gondal`, c: `Marteau`, d: `Yan`, e: `Hu`, f: `Liu`, g: `Hong`, h: `Luoma`, i: `Raufi`) comparing unbounded bulk CLR responder $hat(beta)$ (`Pan-iAtlas-Meta`) against single-cell `milopy` $overline("logFC")$ for the 42 `Consensus-Meta42` states, with OLS regression fits (dashed). The combined normalized overlay of all 9 cohorts is shown in @fig-main-combined[panel c of Figure].
  ]
) <fig-supp-milo-deconv-scatter-grid>

#pagebreak()

#figure(
  image("../figures/milo-deconv/iatlas_classifier_and_reference_tier_pareto.png", width: 100%, height: 178mm, fit: "contain"),
  caption: [
    *Reference-tier Pareto frontier (`#1`–`#13`) and Random Forest non-linear conditional importance outliers.*
    *(a)* Two-objective Pareto landscape plotting cross-modality effect concordance (Pearson $r$ vs. single-cell `milopy` $overline("logFC")$) against multi-cohort predictive `ROC-AUC` across all 13 single-cell references (`#1`–`#13`, colored by Reference Tier 1A, 1B, 2, and 3 and sized by state count $K$), illustrating how `Consensus-Meta42` (`K = 42`) resolves the dimensionality–collinearity penalty of `Unified-All9` (`K = 355`).
    *(b)* Out-of-fold pan-cancer Random Forest Gini feature importance $I_k^("RF")$ regressed against univariate single-cell abundance magnitude $|overline("logFC")_k|$ across all 42 `Consensus-Meta42` states ($N = 1,015$). Six positive-residual outlier states—led by `#1` `M17` (`ISLR+COL5A2+ myCAF`, $I^("RF") = 0.0312$ despite $|overline("logFC")| = 0.45$), `M20` (`IL7R+CD40LG+ Trm`), `M04` (`GNLY+PRF1+ Eff. T`), `M29` (`JCHAIN+MZB1+ Plasma`), `M05` (`TRAF3+CD96+ Mem. T`), and `M30` (`MS4A1+CR2+ TLS B`)—exhibit disproportionate tree-split leverage indicative of non-linear conditional interactions.
  ]
) <fig-pareto-rf-outliers>

#pagebreak()

#figure(
  image("../figures/milo-pairwise/pan_cohort_core_hubs_z_acat_heatmap.png", width: 100%, height: 170mm, fit: "contain"),
  caption: caption_pan_cohort_core_hubs_main
) <fig-pan-core-hubs>

#pagebreak()

#figure(
  image("../figures/milo-pairwise/pan_cohort_unified_pairwise_umap_composite.png", width: 100%, height: 170mm, fit: "contain"),
  caption: caption_pan_cohort_unified_umap_main
) <fig-pan-unified-umap>

#pagebreak()

#figure(
  image("../figures/milo-pairwise/pan_cohort_compartment_pair_umap_grid.png", width: 100%, height: 170mm, fit: "contain"),
  caption: caption_pan_cohort_pair_grid_main
) <fig-pan-pair-grid>

#pagebreak()

// Cross-Modality Deconvolution & Gene-Selection Approaches Comparison Suite
#figure(
  image("../figures/milo-deconv/supp_fig_deconv_factorial_solver_gene_matrix.svg", width: 100%, height: 176mm, fit: "contain"),
  caption: caption_supp_deconv_factorial_matrix
) <fig-supp-deconv-factorial-matrix>

#pagebreak()

#figure(
  image("../figures/milo-deconv/supp_fig_deconv_four_axis_validation_suite.svg", width: 100%, height: 174mm, fit: "contain"),
  caption: caption_supp_deconv_validation_suite
) <fig-supp-deconv-validation-suite>

#pagebreak()

#figure(
  image("../figures/milo-deconv/supp_fig_deconv_approaches_head_to_head_comparison.svg", width: 100%, height: 176mm, fit: "contain"),
  caption: caption_supp_deconv_head_to_head
) <fig-supp-deconv-head-to-head>

#pagebreak()

#figure(
  image("../figures/milo-deconv/supp_fig_deconv_gene_selection_mechanisms.svg", width: 100%, height: 176mm, fit: "contain"),
  caption: caption_supp_deconv_gene_mechanisms
) <fig-supp-deconv-gene-mechanisms>

#pagebreak()

= Extended Supplementary Figures

#set figure(
  supplement: [Extended Supplementary Figure],
  numbering: n => "ES" + str(n),
)
#counter(figure.where(kind: image)).update(0)

// Moved block: formerly S4–S32 (Single-Cell Reference through Unified Factorial Interaction Surface)
#figure(
  image("../figures/single-cell/umap_cell_type_annotations.png", width: 100%, height: 195mm, fit: "contain"),
  caption: [Hierarchical cell type annotation of the single-cell reference. UMAP projection of 40,002 single cells displaying 22 transcriptionally distinct cell types across 6 major lineages.]
) <fig-sc-ref>

#pagebreak()

#figure(
  image("../figures/single-cell/qc_metadata_distributions.png", width: 92%, height: 185mm, fit: "contain"),
  caption: [Single-cell subsampling metadata and QC distributions. Distributions of total unique molecular identifier (UMI) counts, detected genes, mitochondrial percentage, and ribosomal protein percentage across the 40,002 subsampled cells.]
) <fig-supp-s1>

#pagebreak()

#figure(
  image("../figures/single-cell/umap_granularities.png", width: 95%, height: 192mm, fit: "contain"),
  caption: [UMAP projection of Leiden clustering resolutions. Evaluation across clustering resolutions (0.1 to 1.5) on the reference dataset showcasing lineage and subpopulation segregation.]
) <fig-supp-s2>

#pagebreak()

#figure(
  image("../figures/single-cell/umap_tumor_inference.png", width: 95%, height: 192mm, fit: "contain"),
  caption: [UMAP projection of dual-method tumor cell inference. Overlap of epithelial marker gene expression (left) and inferred chromosomal copy number variation (CNV) scores (right) used to identify malignant cells.]
) <fig-supp-s3>

#pagebreak()

#figure(
  image("../figures/single-cell/subclustering_metrics_distributions.png", width: 95%, height: 192mm, fit: "contain"),
  caption: [Daniela Witten Selective Inference p-values. Truncated normal selective inference p-values computed for all sub-clustering configurations, demonstrating statistically significant cluster separation ($p < 10^(-16)$).]
) <fig-supp-s4>

#pagebreak()

// S9–S16: Trial vs. TCGA Baseline Deconvolution, K-Means Classification & RF Prediction
#figure(
  image("../figures/deconvolution/deconv_distributions_kmeans_subcluster_res_0.5.png", width: 100%, height: 195mm, fit: "contain"),
  caption: [Trial versus TCGA baseline deconvolution distributions. Side-by-side horizontal boxplots comparing the 58 sub-cluster fractions between the 9 trial cohorts and matched untreated TCGA indication cohorts.]
) <fig-supp-s5>

#pagebreak()

#figure(
  image("../figures/deconvolution/deconv_split_response_kmeans_subcluster_res_0.5.png", width: 100%, height: 195mm, fit: "contain"),
  caption: [Cell fraction distributions split by clinical response. Boxplot grid displaying K-means 0.5 sub-cluster fractions split by RECIST clinical response (Responders, green vs Non-Responders, red).]
) <fig-supp-s6>

#pagebreak()

#figure(
  image("../figures/deconvolution/deconv_split_tmb_kmeans_subcluster_res_0.5.png", width: 100%, height: 195mm, fit: "contain"),
  caption: [Cell fraction distributions split by tumor mutational burden. Boxplot grid displaying sub-cluster fractions split by cohort median DNA-TMB (High TMB, blue vs Low TMB, orange).]
) <fig-supp-s7>

#pagebreak()

#figure(
  image("../figures/deconvolution/kmeans_cluster_sizes.png", width: 95%, height: 192mm, fit: "contain"),
  caption: [TCGA K-means ($K=100$) cluster sizes and dominant cancer types. Sample count distributions across the 100 K-means centroids colored by primary cancer indication.]
) <fig-supp-s8>

#pagebreak()

#figure(
  image("../figures/deconvolution/kmeans_centroids_heatmap.png", width: 100%, height: 195mm, fit: "contain"),
  caption: [TCGA K-means centroid profiles. Heatmap of the 100 reference centroids across the 58 sub-clusters characterizing primary baseline architectures.]
) <fig-supp-s9>

#pagebreak()

#figure(
  image("../figures/deconvolution/iatlas_classification_umap.png", width: 95%, height: 192mm, fit: "contain"),
  caption: [UMAP embedding of nearest-centroid cancer type classification. Projection of the 1,097 trial samples onto TCGA reference centroids, illustrating the generalization gap in metastatic samples.]
) <fig-supp-s10>

#pagebreak()

#figure(
  image("../figures/deconvolution/deconv_prediction_advanced_grid_kmeans_subcluster_res_0.5.png", width: 100%, height: 195mm, fit: "contain"),
  caption: [Performance metrics (Mean ± SD across 5 random seeds) of Random Forest response prediction using fine-grained deconvolution fractions.]
) <fig-deconv-pred>

#pagebreak()

#figure(
  image("../figures/deconvolution/deconv_prediction_metrics_grid.png", width: 95%, height: 192mm, fit: "contain"),
  caption: [Deconvolution response prediction metrics across feature selection strategies. Comparison of ROC AUC, PR AUC, and MCC across raw features, ANOVA univariate selection, and PCA dimensionality reduction.]
) <fig-supp-s11>

#pagebreak()

// S17: Somatic Mutations & Classifier Overfitting
#figure(
  kind: image,
  grid(
    columns: (1fr, 1fr),
    gutter: 1.2em,
    align(center)[
      #image("../figures/mutations/overfitting_underfitting_rosenberg.png", width: 100%, height: 180mm, fit: "contain")
      #v(0.3em)
      #text(size: 8pt)[(A) Rosenberg Bladder Cohort]
    ],
    align(center)[
      #image("../figures/mutations/overfitting_underfitting_melanoma.png", width: 100%, height: 180mm, fit: "contain")
      #v(0.3em)
      #text(size: 8pt)[(B) Combined Melanoma Cohort]
    ]
  ),
  caption: [Model complexity and cross-validation performance. Comparison of training versus test set ROC AUC across Random Forest, Adaline, and MLP classifiers evaluating genomic overfitting.]
) <fig-supp-s12>

#pagebreak()

// S18–S22: Subclonal VAF, TMB Optimization & Transcriptional Signatures
#figure(
  image("../figures/vaf-tmb/Rosenberg-iAtlas_vaf_tmb_association.png", width: 100%, height: 195mm, fit: "contain"),
  caption: [Subclonal VAF and TMB response association landscape in the Rosenberg bladder cancer cohort, illustrating optimal separation at VAF = 0.05 and TMB = 10.0 mut/Mb.]
) <fig-vaf-peaks>

#pagebreak()

#figure(
  kind: image,
  grid(
    columns: (1fr, 1fr),
    gutter: 1.2em,
    align(center)[
      #image("../figures/vaf-tmb/vaf_density_valley_detection.png", width: 100%, height: 180mm, fit: "contain")
      #v(0.3em)
      #text(size: 8pt)[(A) KDE Valley Detection]
    ],
    align(center)[
      #image("../figures/vaf-tmb/estimated_vs_true_thresholds.png", width: 100%, height: 180mm, fit: "contain")
      #v(0.3em)
      #text(size: 8pt)[(B) Estimated vs True Optimal VAF]
    ]
  ),
  caption: [Unsupervised VAF threshold determination. (A) Kernel Density Estimation valley detection identifying clonal boundaries. (B) Scatter comparison of unsupervised estimated thresholds against empirically optimized cutoffs.]
) <fig-supp-s13>

#pagebreak()

#figure(
  kind: image,
  grid(
    columns: (1fr, 1fr),
    gutter: 1.2em,
    align(center)[
      #image("../figures/vaf-tmb/reliability_score_vs_performance.png", width: 100%, height: 180mm, fit: "contain")
      #v(0.3em)
      #text(size: 8pt)[(A) TRS vs Predictive Strength]
    ],
    align(center)[
      #image("../figures/vaf-tmb/dna_vs_rna_tmb_performance.png", width: 100%, height: 180mm, fit: "contain")
      #v(0.3em)
      #text(size: 8pt)[(B) DNA-TMB vs Expressed RNA-TMB]
    ]
  ),
  caption: [TMB reliability and RNA-TMB performance. (A) Positive correlation (Pearson $r = 0.527$) between TMB Reliability Score (TRS) and predictive strength. (B) Cross-validated ROC AUC comparison between DNA-TMB and Expressed RNA-TMB.]
) <fig-supp-s14>

#pagebreak()

#figure(
  kind: image,
  grid(
    columns: (1fr, 1fr),
    gutter: 1.2em,
    align(center)[
      #image("../figures/vaf-tmb/odds_ratios_forest_plot_MCC.png", width: 100%, height: 180mm, fit: "contain")
      #v(0.3em)
      #text(size: 8pt)[(A) Forest Plot (MCC Objective)]
    ],
    align(center)[
      #image("../figures/vaf-tmb/odds_ratios_forest_plot_F1.png", width: 100%, height: 180mm, fit: "contain")
      #v(0.3em)
      #text(size: 8pt)[(B) Forest Plot (F1 Objective)]
    ]
  ),
  caption: [Forest plots of Odds Ratios for optimized TMB/VAF splits. Fisher's Exact Test Odds Ratios and 95% Confidence Intervals for cutoffs optimized under (A) MCC and (B) F1-Score objectives.]
) <fig-supp-s15>

#pagebreak()

#figure(
  kind: image,
  grid(
    columns: (1fr, 1fr),
    gutter: 1.2em,
    align(center)[
      #image("../figures/signatures/signature_distributions.png", width: 100%, height: 180mm, fit: "contain")
      #v(0.3em)
      #text(size: 8pt)[(A) Signature Score Distributions]
    ],
    align(center)[
      #image("../figures/signatures/signature_overlap_heatmap.png", width: 100%, height: 180mm, fit: "contain")
      #v(0.3em)
      #text(size: 8pt)[(B) Jaccard Overlap Heatmap]
    ]
  ),
  caption: [ImmunoCompass transcriptional signatures. (A) Density distributions of signature scores across trial cohorts. (B) Jaccard similarity heatmap evaluating overlap and redundancy among ImmunoCompass gene sets.]
) <fig-supp-s16>

#pagebreak()

// S23–S26: CD8+ T-Cell, CD4+ T-Cell & Myeloid Single-Cell Landscapes and DGE
#figure(
  image("../figures/immune-subsets/cd8_dge_combined_figure.png", width: 100%, height: 192mm, fit: "contain"),
  caption: [Differential Gene Expression in CD8+ T Cells. (Left) Volcano plot comparing responding versus non-responding patients. (Right) Dot plot of top 15 marker genes characterizing memory persistence versus terminal exhaustion.]
) <fig-cd8-dge>

#pagebreak()

#figure(
  image("../figures/immune-subsets/cd4_dge_combined_figure.png", width: 100%, height: 192mm, fit: "contain"),
  caption: [Differential Gene Expression in CD4+ T Cells. (Left) Volcano plot comparing CD4+ responders versus non-responders. (Right) Dot plot of marker genes illustrating translational machinery vs. DUSP4-mediated senescence.]
) <fig-cd4-dge>

#pagebreak()

#figure(
  kind: image,
  grid(
    columns: (1fr, 1fr),
    gutter: 1.2em,
    image("../figures/immune-subsets/milopy_umap_response_Combined.png", width: 100%, height: 185mm, fit: "contain"),
    image("../figures/immune-subsets/milopy_umap_gradient_Combined.png", width: 100%, height: 185mm, fit: "contain"),
  ),
  caption: [Myeloid spatial landscape. (Left) UMAP colored by clinical response status. (Right) Milo differential abundance log-fold change gradient resolving pro- vs anti-tumor myeloid neighborhoods.]
) <fig-myeloid-spatial>

#pagebreak()

#figure(
  kind: image,
  grid(
    columns: (1fr, 1fr),
    gutter: 1.2em,
    image("../figures/immune-subsets/myeloid_responder_dge_volcano.png", width: 100%, height: 185mm, fit: "contain"),
    image("../figures/immune-subsets/myeloid_responder_dge_dotplot.png", width: 100%, height: 185mm, fit: "contain"),
  ),
  caption: [Differential gene expression in the myeloid compartment. (Left) Volcano plot displaying M1/M2 polarization. (Right) Dot plot of dominant innate marker genes distinguishing responders from non-responders.]
) <fig-myeloid-dge>

#pagebreak()

// S27–S32: Sade-Feldman Matched Clinical Validation & Perturbation Stress-Testing
#figure(
  image("../figures/milo-deconv/sade_feldman_self_concordance.png", width: 100%, height: 188mm, fit: "contain"),
  caption: [Clinical concordance between single-cell Milo DA and pseudobulk deconvolution across 51 melanoma biopsies in the Sade-Feldman cohort. (A) Combined cohort ($n=51$). (B) Pre-treatment baseline cohort ($n=20$). (C) Post-treatment on-therapy cohort ($n=31$). (D) Cohort concordance metrics and ground-truth fidelity. Points colored by diagnostic quadrant: Concordant Responder (green), Concordant Non-Responder (blue), Discordant (red/orange), Neutral (gray).]
) <fig-sade-feldman>

#pagebreak()

#figure(
  image("../figures/milo-deconv/perturbation_sensitivity_curves.png", width: 95%, height: 190mm, fit: "contain"),
  caption: [Multi-replicate single perturbation sensitivity curves ($N=5$ replicates per point). (A) Spearman rank correlation degradation (Mean ± SD). (B) Directional sign concordance (%) with dashed coin-flip baseline. (C) Ground truth deconvolution fidelity ($rho_("true","deconv")$). (D) Maximum degradation ranking sorted by peak correlation loss across 8 single distortion modes.]
) <fig-supp-s17>

#pagebreak()

#figure(
  image("../figures/milo-deconv/synthetic_perturbation_scatter_grid.png", width: 95%, height: 190mm, fit: "contain"),
  caption: [Sade-Feldman style scatter comparison grid for single perturbations. Eight panels displaying Baseline versus 6 distinct single distortion modes at representative moderate strength ($I approx 0.60$), colored by four-quadrant diagnostic classification with linear regression trends and correlation metrics.]
) <fig-supp-s18>

#pagebreak()

#figure(
  image("../figures/milo-deconv/compound_perturbation_sensitivity_curves.png", width: 100%, height: 192mm, fit: "contain"),
  caption: [Multi-factorial compound stress-testing across 5 replicates ($N=5$). (A) Correlation degradation (Mean ± SD). (B) Directional sign concordance. (C) Ground truth fidelity. (D) Clinical severity ranking.]
) <fig-compound-curves>

#pagebreak()

#figure(
  image("../figures/milo-deconv/compound_perturbation_scatter_grid.png", width: 100%, height: 190mm, fit: "contain"),
  caption: [Sade-Feldman style scatter grid comparing Baseline against the 4 compound clinical regimes at representative clinical strength ($I approx 0.6$). (A) Baseline unperturbed. (B) Clinical Core Needle Biopsy. (C) Inflamed TME. (D) SC Technical Noise. (E) Triple Jeopardy. (F) Concordance summary.]
) <fig-compound-scatter>

#pagebreak()

#figure(
  image("../figures/milo-deconv/factorial_2d_interaction_heatmap.png", width: 100%, height: 192mm, fit: "contain"),
  caption: [2D Factorial Interaction Surface between Cell Size Asymmetry and State Activation Confounding ($5 times 5$ grid, $N=5$ replicates per grid cell, 125 total simulations). (A) Mean Spearman rank correlation ($rho$). (B) Mean directional sign agreement (%).]
) <fig-factorial-heatmap>
