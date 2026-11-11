// Auto-generated caption module for Single-Cell Reference Benchmarking Figures
// Figures: Main Summary Figure and Supplementary Figures S52–S56

#let caption_sc9_summary_main = [
  *Multi-cohort single-cell reference integration, clustering resolution, and cell-state granularity govern deconvolution transferability and immunotherapy response prediction.*
  *(a)* Leiden clustering resolution scaling ($"res" = 0.25$ to $2.0$): Out-of-cohort multi-trial objective response predictive capacity (evaluated on $N = 1,015$ patients across 9 solid-tumor ICB cohorts) improves monotonically with fine-grained single-cell clustering resolution (Pearson $r = +0.439, p = 1.2 times 10^(-14)$, Spearman $rho = +0.445, p = 4.7 times 10^(-15)$), confirming that higher resolution separates specialized functional and exhaustion phenotypes obscured at coarse canonical lineage annotations. Red line with overlaid points traces the cross-strategy empirical mean trajectory; shaded boxplots capture distribution medians and interquartile ranges; jittered markers represent distinct reference constructions.
  *(b)* Cell-state granularity scaling: Stratification across seven discrete cell-state regimes ($<15$, $15–19$, $20–24$, $25–29$, $30–34$, $35–39$, and $>= 40$ states) shows steep monotonic scaling with reference dictionary size (Pearson $r = +0.640, p = 1.1 times 10^(-33)$, Spearman $rho = +0.646, p = 1.6 times 10^(-34)$), demonstrating that predictive gain is fundamentally driven by resolving specialized immune effector and stromal barrier sub-states, reaching optimal predictive saturation in the $35–45$ state regime.
  *(c)* Multi-cohort single-cell integration scaling: Number of integrated single-cell reference datasets ($K = 1 dots 9$) regressed against Mean Cohort Response AUC (Pearson $r = +0.251, p = 2.1 times 10^(-5)$, Spearman $rho = +0.229, p = 1.1 times 10^(-4)$), confirming that multi-cohort integration successfully mitigates cohort-specific single-cell batch distortions and yields robust cellular representations that transfer reliably across heterogeneous clinical trials.
  *(d)* Pan-cancer meta-analytic precision–recall benchmark (`Mean PR-AUC` $plus.minus 95%$ CI across $N = 1,015$ pre-treatment standard ICB patients; dashed horizontal line marks the $0.314$ responder prevalence baseline) comparing `Consensus-Meta42` deconvolution inference (`Deconvolution Inference (Our)`, $"PR-AUC" = 0.525$ $[0.465, 0.585]$) against 31 established transcriptomic signatures, single-gene checkpoints, genomic `TMB`, and composite synergy predictors.
  *(e)* Corresponding pan-cancer meta-analytic receiver operating characteristic benchmark (`Mean ROC-AUC` $plus.minus 95%$ CI; dashed horizontal line marks the $0.50$ random baseline) across the same 32 predictors, confirming that `Consensus-Meta42` deconvolution (`Deconvolution Inference (Our)`, $"ROC-AUC" = 0.720$ $[0.669, 0.771]$) significantly outperforms all 31 comparator biomarkers.
]

#let caption_supp_sc9_fig1 = [
  *Exhaustive 18-panel multi-parameter benchmark of single-cell reference predictive capacity across resolution, cell states, multi-dataset integration, and cell subsampling.*
  *(a–c)* Leiden resolution scaling ($"res" = 0.25$ to $2.0$) evaluated on Melanoma ($r = +0.413, p = 6.2 times 10^(-13)$, $rho = +0.407, p = 1.4 times 10^(-12)$), Pan-Cancer pooled strata ($r = +0.380, p = 4.1 times 10^(-11)$, $rho = +0.378, p = 5.2 times 10^(-11)$), and Mean Cohort Response AUC ($r = +0.439, p = 1.2 times 10^(-14)$, $rho = +0.445, p = 4.7 times 10^(-15)$), demonstrating consistent monotonic gains across all evaluation strata.
  *(d–f)* Cell-state granularity regimes ($<15$ to $>= 40$ states) displaying strong positive predictive scaling for Melanoma ($r = +0.573, p = 1.7 times 10^(-25)$, $rho = +0.581, p = 2.7 times 10^(-26)$), Pan-Cancer ($r = +0.596, p = 1.3 times 10^(-28)$, $rho = +0.602, p = 2.8 times 10^(-29)$), and Mean Cohort AUC ($r = +0.640, p = 1.1 times 10^(-33)$, $rho = +0.646, p = 1.6 times 10^(-34)$).
  *(g–i)* Number of integrated single-cell datasets ($K = 1 dots 9$) across Melanoma ($r = +0.187, p = 0.0016$, $rho = +0.163, p = 0.0059$), Pan-Cancer ($r = +0.228, p = 1.1 times 10^(-4)$, $rho = +0.207, p = 4.2 times 10^(-4)$), and Mean Cohort AUC ($r = +0.251, p = 2.1 times 10^(-5)$, $rho = +0.229, p = 1.1 times 10^(-4)$), demonstrating positive transfer scaling with dataset integration depth.
  *(j–l)* Continuous log-cell scaling curves ($N_("cells") = 1,000$ to $> 200,000$ cells) confirming logarithmic performance saturation with increasing cell count (Melanoma $r_("log") = +0.401$, Pan-Cancer $r_("log") = +0.418$, Mean Cohort $r_("log") = +0.432$, all $p < 10^(-10)$).
  *(m–o)* Cell subsampling titration ($10%$ to $100%$ downsampling) displaying performance stability down to $25%$ subsampling before experiencing modest attenuation.
  *(p–r)* Reference strategy classification (Single Dataset vs. Criteria-Combined vs. Random-Combined) and marginal AUC gains per resolution step ($Delta "AUC"$), showing positive marginal returns up to $"res" = 1.50$ followed by asymptotic stabilization.
]

#let caption_supp_sc9_fig2 = [
  *Deconvolution signature collinearity, condition number scaling, variance inflation, and shrinkage regularization across single-cell reference granularities.*
  *(a–c)* Signature matrix condition number $kappa(bold(X))$ regressed against Leiden clustering resolution across single-cell reference cohorts, revealing exponential condition number growth ($kappa > 10^3$) at high resolutions ($"res" >= 1.50$) driven by shared transcriptomic programs among closely related lymphocyte and myeloid sub-lineages.
  *(d–f)* Maximum Variance Inflation Factor ($"VIF"_("max")$) across reference cell states as a function of dictionary size $K$, isolating severe multi-collinearity in naive unregularized regression and highlighting the necessity of $ell_2$ shrinkage in `InstaPrism` and Bayesian prior constraints in `BayesPrism`.
  *(g–i)* Pairwise cell-state signature Pearson correlation distributions and condition number mitigation via automated thresholding and ridge shrinkage, establishing the mathematical bounds of stable deconvolution in high-dimensional microenvironments.
]

#let caption_supp_sc9_fig3 = [
  *Multi-criteria reference selection decision framework, Pareto trade-off frontiers, and optimal operating regimes.*
  *(a–c)* Two-dimensional Pareto trade-off surfaces plotting cross-cohort predictive accuracy (Mean Cohort AUC) against signature condition number penalty ($log_10(kappa)$), delineating the non-dominated Pareto frontier across 40+ evaluated reference constructions.
  *(d–f)* Harmonic decision score and multi-objective optimization balancing predictive power, biological interpretability, computational complexity, and collinearity risk, identifying the $K = 35–45$ state regime as the global optimum for clinical deconvolution.
  *(g–h)* Sensitivity analysis and decision framework boundary mapping across individual tumor histologies (melanoma, bladder, renal, and gastrointestinal cancers).
]

#let caption_supp_sc9_fig4 = [
  *Winning classifier architecture, biomarker driver weights, and cross-cohort biological stability across high-resolution reference states.*
  *(a–c)* Cross-validated predictive performance across machine learning paradigms (`Logistic L2`, `Linear Reg.`, `Random Forest`, `Elastic-Net`, `Lasso L1`) on the winning reference configurations, highlighting regularized linear models for small cohorts and non-linear tree ensembles for pooled pan-cancer cohorts.
  *(d–f)* Normalized model coefficients and feature importance rankings identifying universal responder-associated cellular states (anchored by `CXCL13+PDCD1+` exhausted CD8+ T cells and `CXCL9+STAT1+` M1-TAMs) and non-responder-associated barrier states (anchored by `FAP+` mature CAFs, `S100A8+` MDSCs, and mesenchymal EMT tumor cells).
  *(g–i)* Cross-cohort biological stability and sign concordance across independent clinical trial validations, confirming consistent mechanistic polarity across human solid tumors.
]

#let caption_supp_sc9_fig5 = [
  *Selected individual cohort AUC response profiles, marginal integration gains, and statistical slope significance tests.*
  *(a–d)* Individual clinical trial predictive response profiles across single-cell reference integration depths for representative cohorts: Melanoma (`Gide`, `Hugo`, `Riaz`), Bladder (`Rosenberg`), Renal Cell Carcinoma (`McDermott`, `Choueiri`), and Pancreatic/TNBC (`Padron`), displaying positive response trends across distinct histological contexts.
  *(e–h)* Marginal predictive gain ($Delta "ROC-AUC"$) per integrated reference dataset alongside empirical linear slope estimates with exact Wald $p$-values and within-cohort permutation tests, confirming statistically significant positive performance slopes ($p < 0.05$ marked by asterisks) across all major lineages.
  *(i)* Pan-cancer composite summary of marginal integration returns, establishing that cross-cohort reference harmonization yields statistically verifiable improvements in out-of-cohort clinical response prediction.
]
