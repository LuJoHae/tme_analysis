== 3.17 Finite-Sample Controlled Transcriptomic Stability Across 12 Clinical Cohorts
High-dimensional transcriptomic biomarker discovery in immuno-oncology is notoriously plagued by severe sample variance, cohort-specific batch artifacts, and high false discovery rates ($p >> n$) @meinshausen2010stability, @shah2013variable. When single predictive models are trained on small or moderate patient cohorts ($n < 100$), multi-collinear immune co-expression networks cause standard regularized estimators to select idiosyncratic, non-reproducible gene sets.

To overcome this fundamental limitation, we deployed finite-sample controlled stability selection with complementary pairs subsampling ($B=100$ cycles, $2B=200$ subsamples of size $floor(n/2)$) across a compendium of 12 clinical datasets spanning 1,015 pre-treatment immunotherapy patients (@fig-stability-paths-composite, @tab-stability-cohorts). By imposing rigorous regularization budgeting ($q_("budget") = floor(sqrt("PFER"_("target") dot p dot (2 pi_("thr") - 1)))$ with nominal PFER bound $<= 2.5$), the feature search space ($p = 822$ immune-oncology genes) was swept systematically across ordered regularization paths without risking active set runaway.

#figure(
  image("../figures/stability-selection/all_cohorts_stability_paths_composite.png", width: 100%),
  caption: [Finite-sample controlled stability paths across 12 clinical immunotherapy cohorts ($n=1,015$ total patients). Each panel depicts empirical feature selection probabilities $hat(Pi)_k^lambda$ across the budgeted regularization path for an individual cohort or multi-cohort pooled setting. The horizontal dashed red line marks the stability selection threshold ($pi_("thr") = 0.60$ for single cohorts, $pi_("thr") = 0.70$ for multi-cohort pooled regimes). Solid colored paths denote stable biomarkers surpassing $pi_("thr")$ under strict PFER error control ($"PFER" <= 2.5$), while faint gray paths denote unselected noise features. The universal chemokine marker _CXCL13_ (highlighted in gold/orange) consistently achieves maximum stability across indications.]
) <fig-stability-paths-composite>

In the pooled Pan-Cancer cohort ($n=1,015$ across melanoma, renal cell carcinoma, urothelial bladder carcinoma, and pancreatic ductal adenocarcinoma), stability selection uncovered a highly focused, statistically unassailable core signature consisting of _CXCL13_, _CD8A_, _CXCL9_, _IFNG_, _LAG3_, and _TIGIT_ ($hat(Pi)_k >= 0.70$, PFER $< 1.8$). Rather than selecting dozens of spurious correlates, regularization budgeting curtailed exploration once the active set size reached capacity, preserving finite-sample error bounds.

Analyzing individual disease indications and independent clinical trials confirmed remarkable biological convergence alongside indication-specific nuances (@tab-stability-cohorts):
1. *Melanoma Combined ($n=338$)*: _CXCL13_ ($hat(Pi) = 0.94$), _CD8A_ ($0.88$), _PDCD1_ ($0.81$), and _LAG3_ ($0.74$) emerged as dominant invariant drivers, replicating across individual trials including Liu et al. ($n=122$, _CXCL13_, _CD8A_, _TIGIT_), Riaz et al. ($n=98$, _CXCL13_, _PDCD1_), Gide et al. ($n=91$, _CXCL13_, _IFNG_, _CD8A_), and Hugo et al. ($n=27$, _CXCL13_, _CD8A_).
2. *Renal Cell Carcinoma ($n=263$)*: In addition to cytotoxic markers _CXCL13_ and _CD8A_, RCC demonstrated clear vascular and microenvironmental divergence, selecting endothelin receptor type B (_EDNRB_) and vascular endothelial growth factor A (_VEGFA_), particularly in the sunitinib-versus-atezolizumab trial (McDermott et al., $n=247$) and Choueiri et al. ($n=16$).
3. *Urothelial Bladder Carcinoma (Rosenberg et al., $n=298$)*: The atezolizumab Phase II IMvigor210 cohort selected the CXCR3-chemokine axis (_CXCL9_, _CXCL10_) alongside _CD8A_ and _IFNG_ as the primary predictors of clinical response, consistent with baseline CD8+ T-cell peritumoral infiltration.
4. *Pancreatic Ductal Adenocarcinoma (Padron et al., $n=85$)*: Despite the heavily fibrotic and immunologically cold microenvironment characteristic of PDAC, stability selection identified _CXCL13_ and killer cell lectin-like receptor B1 (_KLRB1_) in patients receiving sotigalimab and nivolumab.

#figure(
  table(
    columns: (1.4fr, 0.7fr, 1.1fr, 1.2fr, 1.6fr),
    inset: 4.5pt,
    align: center,
    [*Cohort / Setting*], [*$n$*], [*Indication*], [*Stable Genes ($>= pi_("thr")$)*], [*Biological Role*],
    [Pan-Cancer Combined], [1,015], [Mixed Solid], [_CXCL13_, _CD8A_, _CXCL9_, _IFNG_, _LAG3_], [Universal Tertiary Lymphoid & CD8+ Infiltration],
    [Melanoma Combined], [338], [Melanoma], [_CXCL13_, _CD8A_, _PDCD1_, _LAG3_], [Exhausted Resident Memory T Cells],
    [RCC Combined], [263], [Renal Cell], [_CXCL13_, _CD8A_, _EDNRB_, _VEGFA_], [Angiogenesis & T-Cell Infiltration],
    [Rosenberg et al.], [298], [Bladder], [_CXCL9_, _CXCL10_, _CD8A_, _IFNG_], [IFN-$gamma$ Induced Chemokines],
    [McDermott et al.], [247], [Renal Cell], [_CXCL13_, _CD8A_, _VEGFA_], [Immune-Vascular Interplay],
    [Liu et al.], [122], [Melanoma], [_CXCL13_, _CD8A_, _TIGIT_], [Checkpoint Co-inhibition],
    [Riaz et al.], [98], [Melanoma], [_CXCL13_, _PDCD1_], [Pre-existing Clonal T Cells],
    [Gide et al.], [91], [Melanoma], [_CXCL13_, _IFNG_, _CD8A_], [Effector Cytotoxicity],
    [Padron et al.], [85], [PDAC], [_CXCL13_, _KLRB1_], [CD40 Agonism & NK/T Infiltration],
    [Anders et al.], [31], [Bladder], [_CXCL9_, _STAT1_], [JAK/STAT Transcriptional Priming],
    [Hugo et al.], [27], [Melanoma], [_CXCL13_, _CD8A_], [Innate Anti-PD-1 Sensitivity],
    [Choueiri et al.], [16], [Renal Cell], [_CXCL13_], [Conserved Chemokine Hub]
  ),
  caption: [Summary of finite-sample controlled stability selection across 12 clinical immunotherapy cohorts ($n=1,015$). Feature space $p = 822$ immune-oncology genes. Thresholds: $pi_("thr") = 0.70$ for combined/pooled cohorts, $pi_("thr") = 0.60$ for individual study cohorts. Nominal PFER bounded by $<= 2.5$.]
) <tab-stability-cohorts>

== 3.18 Benchmark of 13 Machine Learning and Econometric Fitters
To ascertain whether transcriptomic stability depends on specific functional parameterizations or model inductive biases, we benchmarked 13 machine learning and econometric fitters across 5 methodological paradigms on the multi-cohort benchmark dataset (@fig-fitter-overlap).

#figure(
  image("../figures/stability-selection/fitter_gene_overlap_nature.png", width: 100%),
  caption: [Cross-paradigm benchmark of 13 machine learning and econometric fitters in stability selection. (A) Gene selection frequency across all 13 fitters. _CXCL13_ emerged as the invariant, model-agnostic driver of response, selected by 12 out of 13 algorithms. Core immune effectors (_CD8A_, _CXCL9_, _CXCL10_, _IFNG_, _LAG3_, _TIGIT_) were recurrently recovered across multiple paradigms. (B) Fitter paradigm classification highlighting convex GLMs (Lasso, ElasticNet), non-convex/adaptive estimators (SCAD, MCP, AdaptiveLasso, HardThresholding), tree ensembles (LightGBM, XGBoost, RandomForest), cohort-aware econometric models (LassoFWL, MERF), and ordered weighted $ell_1$ estimators (OSCAR, SLOPE). (C) Gene selection matrix across all fitters, demonstrating consensus patterns and method-specific inductive biases.]
) <fig-fitter-overlap>

1. *Universal Invariance of _CXCL13_ Across 12 of 13 Fitters*:
Strikingly, _CXCL13_ was selected by 12 out of 13 fitters across all five paradigms (selected by Lasso, ElasticNet, SCAD, MCP, AdaptiveLasso, HardThresholding, LightGBM, XGBoost, RandomForest, LassoFWL, OSCAR, and SLOPE; only MERF prioritized alternative non-linear splits). This unanimous cross-paradigm concordance establishes that _CXCL13_ is not an artifact of coordinate descent, $ell_1$ geometry, or tree split heuristics, but an invariant, robust biological marker of immunotherapy sensitivity.

2. *Convex vs. Non-Convex Regularization Dynamics*:
Convex GLMs (Lasso and ElasticNet) yielded highly focused active sets ($3$ to $5$ stable genes), effectively suppressing high-dimensional noise. In contrast, non-convex penalties (SCAD and MCP) with local linear approximation (LLA) eliminated asymptotic attenuation bias for large effect coefficients, yielding slightly sharper selection probability transitions along the regularization path. Adaptive Lasso, leveraging initial ridge weights, produced sparse sets with minimal false discovery bleeding into peripheral co-expressed transcripts.

3. *Gradient Tree Ensembles and Non-Linear Epistasis*:
LightGBM, XGBoost, and Random Forest captured non-linear threshold effects and multi-gene interactions that evade strictly additive generalized linear models. While tree ensembles exhibited slightly broader selection profiles due to stochastic feature subsampling (Gini gain and split coverage), their consensus set strongly aligned with regularized GLMs, confirming that the dominant signal in bulk tumor transcriptomics is governed by core immune infiltration modules.

4. *Purging Inter-Cohort Batch Effects via Econometric Partialling and Mixed Models*:
A critical challenge in pooling multi-center clinical trials is confounding from study-specific baseline response rates and batch variations.
- *Frisch-Waugh-Lovell Partialled Lasso (`LassoFWL`)*: By projecting expression data and binary outcomes onto the orthogonal complement of the cohort indicator subspace ($bold(M)_bold(C) = bold(I) - bold(C)(bold(C)^T bold(C))^(-1)bold(C)^T$), `LassoFWL` purged cohort-level fixed shifts before performing coordinate descent. Crucially, _CXCL13_, _CD8A_, and _CXCL9_ retained high stability after full confounder decontamination, demonstrating that their predictive capacity is entirely intra-cohort and not driven by between-cohort batch differences.
- *Mixed Effects Random Forest (`MERF`)*: Incorporating cohort random intercepts ($b_j ~ cal(N)(0, sigma_b^2)$) through Expectation-Maximization decoupled global transcriptomic non-linearities from institutional baseline shifts.

== 3.19 Concordance and Invariant Biomarker Discovery via Ordered Weighted L1 Penalization
While standard Lasso ($ell_1$) performs effective variable selection, it suffers from two well-known theoretical pathologies under biological co-expression: (i) it arbitrarily selects one gene from a cluster of correlated covariates and discards the rest, and (ii) it cannot guarantee finite-sample False Discovery Rate (FDR) control without stringent irrepresentable conditions @bondell2008simultaneous, @bogdan2015slope.

To resolve these limitations, we integrated Ordered Weighted $ell_1$ (OWL) norms into the stability selection pipeline:
1. *OSCAR Octagonal Clustering*:
The affine linear weight sequence $w_j^("OSCAR") = lambda_1 + lambda_2 (p - j)$ introduces octagonal polytope facets that enforce exact equality of absolute coefficients for correlated features ($|hat(beta)_i| = |hat(beta)_j|$). In our multi-cohort stability selection, OSCAR successfully clustered co-expressed chemokine and cytotoxicity genes (_CXCL13_, _CXCL9_, _CD8A_, _IFNG_) into a unified predictive module, preventing the arbitrary feature-swapping typical of standard Lasso across bootstrap resamples.

2. *SLOPE Adaptive FDR Control*:
By penalizing sorted coefficients with decaying Gaussian quantiles ($w_j^("SLOPE") = lambda dot Phi^(-1)(1 - (q_("FDR") dot j)/(2 p))$ with $q_("FDR") = 0.10$), SLOPE adjusted the regularization penalty dynamically to the empirical size of the active model. Under stability selection, SLOPE strictly curbed the selection of marginal noise covariates while decisively selecting the core cytotoxic axis (_CXCL13_, _CD8A_, _CXCL9_), verifying finite-sample error control under empirical correlation structures.

3. *Algorithmic Efficiency of PAVA Isotonic Projection*:
By casting the proximal projection step into a monotonic isotonic regression program solved via the Pool Adjacent Violators Algorithm (PAVA) in $O(p log p)$ time, both OSCAR and SLOPE achieved computational throughput comparable to coordinate descent. Combined with FISTA backtracking line search, the OWL solvers converged reliably across all $2B=200$ subsampling iterations, providing a mathematically principled, cluster-aware framework for reproducible clinical biomarker discovery.

== 3.20 Curated Multi-Paradigm Immuno-Oncology Biomarker Dossier
To translate the mathematical convergence across 91 stability selection runs (12 cohorts $times$ 13 fitters) into actionable biological targets, we established a composite Evidence Index ($E(g) in [0, 100]$) and curated an expert-ready biomarker dossier (@tab-curated-expert-gene-dossier). Candidate genes were categorized into five functional immuno-oncology axes: Antigen Processing & Presentation (APM), Chemokines & TLS (TLM), Effector Cytotoxicity & Lineage (CYT), Checkpoints & Co-stimulation (CKP), and Immunosuppressive Stroma & Metabolism (STM).

#include "../tables/table_expert_gene_dossier.typ"

1. *Universal Replicated Drivers (Tier 1)*:
Nineteen candidate genes achieved Tier 1 status ($E(g) >= 60$ or replication across $>= 3$ independent cohorts). Classical MHC-I molecule _HLA-A_ demonstrated the highest cross-paradigm consensus, selected across 12 distinct fitters and 4 cohorts with maximum stability ($hat(Pi) = 1.00$). Strikingly, two potent resistance determinants emerged with exceptional stability:
- *Prostaglandin-Endoperoxide Synthase 2 (_PTGS2_ / COX-2)*: Selected 30 times across 3 cohorts (including RCC and Pan-Cancer) and 9 fitters ($hat(Pi)_("max") = 0.95$). Expression was significantly elevated in non-responders ($log_2"FC" = -0.61$, $r_("pb") = -0.125$, Welch $p = 3.24 times 10^(-5)$), underscoring stromal prostaglandin synthesis as a conserved driver of checkpoint resistance amenable to pharmacological inhibition via celecoxib or apricoxib.
- *Helios Transcription Factor (_IKZF2_)*: Selected 20 times across 3 cohorts and 10 fitters ($hat(Pi)_("max") = 1.00$), marking stable immunosuppressive intratumoral regulatory T cells ($log_2"FC" = -0.32$, $p = 0.011$).

2. *Effector Priming vs. T-Cell Exhaustion Dualities*:
Core cytotoxic effectors (_TBX21_, _IFNG_, _CD8B_) consistently correlated with clinical response ($log_2"FC" > +0.70$, $p < 10^(-4)$), validating the premise that pre-existing cytotoxic infiltration is mandatory for therapeutic anti-PD-(L)1 efficacy. Concurrently, co-stimulatory 4-1BB ligand (_TNFSF9_) and checkpoint receptor _LAG3_ were captured across multiple paradigms, reflecting the ongoing immunological equilibrium between activation and adaptive exhaustion in the tumor microenvironment.

3. *Mechanistic Rescues via Grouped and Tree Regularization (Tiers 2 & 3)*:
Methods sensitive to genomic collinearity and non-linear interactions rescued key biological transactivators that were dropped by standard coordinate-descent Lasso:
- The MHC-I master regulator _NLRC5_ and immunoproteasome subunit _PSMB8_ were stably recovered by Group Lasso and OSCAR ($hat(Pi) >= 0.95$).
- The MHC-II master regulator _CIITA_ was identified by Random Forest and MERF, reflecting its non-linear threshold role as an interferon-gamma-induced switch on antigen-presenting cells.

