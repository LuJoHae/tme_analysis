== 3.10 Single-Cell Neighborhood Abundance versus Pseudobulk Deconvolution
A central unresolved question in computational oncology is whether bulk transcriptomic deconvolution faithfully recovers the specific cell states identified by single-cell differential abundance (DA) testing @dann2022differential, @chu2022cell. While both methodologies aim to identify microenvironmental subsets predictive of clinical outcome, they operate upon fundamentally distinct mathematical representations:

1. *Count Abundance on Graph Manifolds vs. mRNA Mass Fractions*: Milo tests physical cell numbers $N_(k,j)$ within localized, overlapping kNN graph neighborhoods using Negative Binomial generalized linear models. In contrast, bulk deconvolution algorithms (such as `instaprism` and `BayesPrism`) invert linear mixtures on the simplex $Delta^(K-1)$, inherently estimating mRNA mass fractions $f_(k,j)^("mRNA")$. When cell types exhibit unequal transcriptional outputs (e.g. metabolically active macrophages or plasma cells producing $10times - 50times$ more total mRNA per cell than resting lymphocytes, $S_k >> S_m$), the mass fraction decouples from physical cell abundance:
$ f_(k,j)^("mRNA") = (p_(k,j) S_k) / (sum_m p_(m,j) S_m) $
2. *Local Graph Geometry vs. Ill-Conditioned Matrix Inversion*: Milo operates non-parametrically across continuous phenotypic gradients without requiring discrete cluster boundaries. Deconvolution, however, requires a discrete signature matrix $bold(Phi) in bb(R)^(G times K)$. When closely related sub-lineages share substantial marker programs (high reference collinearity), the condition number $kappa(bold(Phi))$ inflates, causing matrix inversion instability and sign-flipping between collinear cell states.
3. *Abundance vs. Transcriptional Activation Confounding*: If an immune cell state becomes activated in responding patients (upregulating effector transcripts $5times - 20times$ without cellular proliferation), Milo detects constant cell density on the kNN graph ($overline("logFC") approx 0$). Deconvolution, however, detects a surge in total bulk marker reads, falsely inferring massive cellular expansion ($hat(beta)^("deconv") >> 0$).

== 3.11 Clinical Concordance in the Sade-Feldman Melanoma Cohort
To establish an empirical clinical benchmark, we evaluated matched single-cell Milo DA and pseudobulk deconvolution across 51 pre- and post-treatment tumor biopsies from the Sade-Feldman metastatic melanoma cohort (GSE120575, 12 immune states) @sadefeldman2018defining (@fig-sade-feldman).

#figure(
  image("../figures/milo-deconv/sade_feldman_self_concordance.png", width: 95%),
  caption: [Clinical concordance between single-cell Milo DA and pseudobulk deconvolution across 51 melanoma biopsies in the Sade-Feldman cohort. (A) Combined cohort ($n=51$). (B) Pre-treatment baseline cohort ($n=20$). (C) Post-treatment on-therapy cohort ($n=31$). (D) Cohort concordance metrics and ground-truth fidelity. Points colored by diagnostic quadrant: Concordant Responder (green), Concordant Non-Responder (blue), Discordant (red/orange), Neutral (gray).]
) <fig-sade-feldman>

In the combined cohort, Milo effect sizes and deconvolution regression coefficients exhibited strong, highly significant rank correlation (*Spearman $rho = 0.853$*, $p = 4.18 times 10^(-4)$; Pearson $r = 0.837$, $p = 6.92 times 10^(-4)$), with *75.0%* directional sign concordance (10 of 12 states agreeing). High concordance persisted in pre-treatment baseline biopsies ($rho = 0.818$, $p = 1.14 times 10^(-3)$, 66.7% sign agreement) and post-treatment biopsies ($rho = 0.783$, $p = 2.59 times 10^(-3)$, 75.0% sign agreement).

Classifying cell states into diagnostic quadrants identified consistent clinical biomarkers (@tab-sade-feldman):
- *Concordant Responders*: `04_Naive B cells` emerged as the single most potent concordant responder biomarker (Milo $overline("logFC") = +3.635$, deconvolution $hat(beta) = +2.066$, Wald $p = 3.6 times 10^(-4)$), corroborating the critical role of tertiary lymphoid structures in melanoma immunotherapy response.
- *Concordant Non-Responders*: `05_Macrophages` ($hat(beta) = -1.562$, $p = 0.0039$), `07_Tem/Trm cytotoxic T cells` ($hat(beta) = -1.047$, $p = 0.0395$), and `10_pDC` ($hat(beta) = -0.237$) consistently mapped to the negative response quadrant.
- *Collinear Discordance*: Strikingly, `01_Tem/Trm cytotoxic T cells` displayed marked discordance in pre-treatment biopsies (Milo $overline("logFC") = +1.726$, deconvolution $hat(beta) = -0.509$). Because 5 of the 12 clusters represented closely related cytotoxic CD8+ subsets, flat matrix inversion subtracted fraction from `01_Tem/Trm` to mathematically offset the depletion of collinear `07_Tem/Trm` due to negative covariance ($kappa(bold(Phi)^("state")) = 1,248.6$). Applying BayesPrism's hierarchical aggregation operator $bold(M)$ to consolidate these 5 sibling subsets into a unified `Cytotoxic_T` lineage annihilated the collinear noise vector ($bold(v) in "Null"(bold(M))$), reducing the condition number by 68.2% ($kappa = 397.1$) and restoring positive concordance ($hat(beta) = +0.814$).

#figure(
  table(
    columns: (1.5fr, 1fr, 1.1fr, 1.1fr, 1.6fr),
    inset: 4.5pt,
    align: center,
    [*Cell State*], [*Milo logFC*], [*Deconv $hat(beta)$*], [*Deconv $p$*], [*Diagnostic Quadrant*],
    [04_Naive B cells], [+3.635], [+2.066], [0.00036], [Concordant Responder],
    [06_Tem/Temra], [-0.120], [+0.173], [0.718], [Discordant (Milo-, Deconv+)],
    [08_Treg], [+2.402], [+0.068], [0.887], [Concordant Neutral],
    [10_pDC], [-1.058], [-0.237], [0.621], [Concordant Non-Responder],
    [01_Tem/Trm], [+1.726], [-0.509], [0.294], [Discordant (Milo+, Deconv-)],
    [00_Tem/Trm], [-0.241], [-0.520], [0.284], [Concordant Non-Responder],
    [09_Plasma cells], [-2.318], [-0.743], [0.133], [Concordant Non-Responder],
    [03_Tem/Trm], [-1.348], [-0.820], [0.099], [Concordant Non-Responder],
    [07_Tem/Trm], [-1.324], [-1.047], [0.0395], [Concordant Non-Responder],
    [05_Macrophages], [-1.365], [-1.562], [0.0039], [Concordant Non-Responder]
  ),
  caption: [Biomarker state classification in pre-treatment melanoma biopsies (Sade-Feldman GSE120575).]
) <tab-sade-feldman>

== 3.12 Multi-Replicate Perturbation Benchmarking Identifies Primary Failure Modes
To map the exact failure thresholds of deconvolution, we conducted an automated multi-replicate benchmark evaluating 8 distinct perturbation modes across 5 random seeds per condition ($N=5$, 240 simulations, @fig-supp-s17, @fig-supp-s18).

This systematic stress-testing uncovered two catastrophic single-variable failure modes:
1. *Patient-Specific Marker Dysregulation is the Single Most Lethal Biological Distortion*: When tumors suppress or downregulate lineage markers in bulk tissue (e.g. HLA downregulation, hypoxia-induced exhaustion pruning; $sigma_("dysreg") in [0, 2.0]$), deconvolution Spearman $rho$ plummeted from *$+0.967$ down to $-0.043 plus.minus 0.473$* at extreme intensity, while directional sign concordance collapsed below a coin-flip (*$47.5 plus.minus 18.5%$*), and ground-truth fidelity collapsed to *$-0.062 plus.minus 0.431$*. Because the deconvolution reference matrix $bold(Phi)$ is static, marker suppression in a patient's bulk tissue is interpreted as complete cell absence, entirely decoupling from physical single-cell abundance.
2. *Cell Size Asymmetry Drives Severe Rank Degradation*: Scaling per-cell mRNA content ($S_k in [1times, 50times]$) collapsed Spearman $rho$ from *$+0.824$ to $+0.252 plus.minus 0.151$* ($Delta rho = 0.572$), with deconvolution fidelity dropping to *$+0.257 plus.minus 0.072$*. High-mRNA lineages (such as macrophages and plasma cells) systematically monopolized bulk reads, distorting lymphocyte proportions.
3. *State Activation Confounding Inverts Directional Signs*: Elevated marker transcription ($gamma in [1times, 20times]$) dropped directional sign agreement from *82.5%* down to *42.5%* ($I=0.4$), demonstrating that cytokine-driven gene induction can invert perceived response associations in bulk deconvolution.

== 3.13 Multi-Factorial Compound Regimes & 2D Interaction Surface
Because clinical biopsies experience simultaneous distortions, we benchmarked 4 compound clinical regimes ($N=5$ replicates, 120 simulations, @fig-compound-curves, @fig-compound-scatter):

#figure(
  image("../figures/milo-deconv/compound_perturbation_sensitivity_curves.png", width: 95%),
  caption: [Multi-factorial compound stress-testing across 5 replicates ($N=5$). (A) Correlation degradation (Mean ± SD). (B) Directional sign concordance. (C) Ground truth fidelity. (D) Clinical severity ranking.]
) <fig-compound-curves>

#figure(
  image("../figures/milo-deconv/compound_perturbation_scatter_grid.png", width: 95%),
  caption: [Sade-Feldman style scatter grid comparing Baseline against the 4 compound clinical regimes at representative clinical strength ($I approx 0.6$). (A) Baseline unperturbed. (B) Clinical Core Needle Biopsy. (C) Inflamed TME. (D) SC Technical Noise. (E) Triple Jeopardy. (F) Concordance summary.]
) <fig-compound-scatter>

The *Clinical Core Needle Biopsy* regime (Cell Size + Sampling Sparsity + Ghost Contamination) emerged as the single most destructive condition across all compound benchmarks:
- Spearman $rho$ degraded from *$+0.943 plus.minus 0.055$* down to *$+0.271 plus.minus 0.265$*.
- Directional sign concordance collapsed to *52.5%* (pure coin-flip).
- Deconvolution ground-truth fidelity crashed to *$+0.229 plus.minus 0.319$*.

Conversely, regimes driven primarily by collinearity and patient-level expression shifts (*Inflamed TME*, $rho = +0.967$; *SC Technical Noise*, $rho = +0.776$) retained high rank stability because Milo's generalized linear model successfully controlled for patient batch variance, and `instaprism`'s Bayesian shrinkage prior buffered against ill-conditioned matrix inversion.

#figure(
  image("../figures/milo-deconv/factorial_2d_interaction_heatmap.png", width: 95%),
  caption: [2D Factorial Interaction Surface between Cell Size Asymmetry and State Activation Confounding ($5 times 5$ grid, $N=5$ replicates per grid cell, 125 total simulations). (A) Mean Spearman rank correlation ($rho$). (B) Mean directional sign agreement (%).]
) <fig-factorial-heatmap>

Finally, evaluating the 2D Factorial Interaction Surface ($5 times 5$ grid of Cell Size $times$ Activation, $N=125$ simulations, @fig-factorial-heatmap) revealed fundamental orthogonal failure mechanics:
- *Cell Size Asymmetry* is the dedicated driver of *rank correlation breakdown*, driving $rho$ from $+0.96$ down to $+0.46$ along the horizontal axis.
- *State Activation Confounding* is the dedicated driver of *directional sign discordance*, causing clusters to drop from $>95%$ down toward $65% - 75%$ along the vertical axis.
- Together, these findings demonstrate that bulk deconvolution and single-cell graph DA are concordant in high-purity, uniform tissue, but decouple under low biopsy yield, unmodeled parenchymal tumor, and microenvironmental marker dysregulation.

== 3.14 Resolving Collinearity Breakdown via Graph-Laplacian Regularization on the Probability Simplex
Having demonstrated that broad lineage aggregation (Section 2.8) stabilizes deconvolution by projecting out collinear variance directions ($bold(v) in "Null"(bold(M))$), we next asked whether fine-grained cell states can be resolved directly on the simplex without collapsing states into coarse types. To systematically evaluate this, we established a controlled continuous-spectrum synthetic benchmark spanning from completely orthogonal states ($r_("sibling") approx 0.0$) through moderate, high, and severe collinearity up to near-singular profiles ($r_("sibling") = 0.990$) across $S=12$ purely synthetic in silico states grouped into $T=4$ lineages ($N=50$ mixture samples per grid point, depth $= 80,000$ reads) (@fig-synthetic-datasets, @fig-collinearity-benchmark, @fig-ground-truth-comparison, @tab-collinearity-benchmark). All reference signatures, proportions, and mixtures are generated from first principles via Gamma and Dirichlet distributions without drawing from empirical single-cell counts.

#figure(
  image("../figures/deconvolution/synthetic_dataset_overview.png", width: 95%),
  caption: [Controlled synthetic benchmark dataset architecture. (A) Reference signature correlation matrices across sibling correlation levels ($r = 0.0, 0.6, 0.9, 0.99$) across $S=12$ synthetic states and $T=4$ lineages (`Lineage_1` to `Lineage_4`), demonstrating emerging $3 times 3$ intra-lineage block collinearity while preserving cross-lineage orthogonality ($r approx 0.0$). (B) Ground truth non-uniform Dirichlet proportion distribution across the 12 synthetic subsets (`L1_State_1` through `L4_State_3`, $N=50$ samples, $bold(alpha) = [1.2, 1.0, 0.8, 1.0, 0.8, 0.6, 1.2, 1.0, 0.8, 0.8, 0.6, 0.4]$). (C) Representative bulk mixture sequencing read counts (Sample 01, depth $= 80,000$ reads across synthetic genes).]
) <fig-synthetic-datasets>

#figure(
  image("../figures/deconvolution/collinearity_regularized_synthetic_benchmark.png", width: 95%),
  caption: [Continuous-spectrum statistical benchmark of collinearity-aware regularized deconvolution across sibling correlations ($r in [0.0, 0.99]$, $S=12$ states, $N=50$ samples per point) across all seven methods: unregularized NNLS, RegDeconv, Rectangle (DWLS-QP), CIBERSORT (reimplementation, $nu$-SVR), InstaPrism, BayesPrism (Gibbs), and the official CIBERSORTx container (`cibersortx/fractions:latest`). (A) Reference Hessian condition number $kappa(bold(H))$ on $log_(10)$ scale: standard Poisson likelihood escalates toward $3.8 times 10^3$ as $r -> 0.99$, whereas Adaptive Deficit Regularization bounds $kappa <= 2,000$ ($kappa = 1,618$). (B) State-level proportion Mean Squared Error (MSE, Mean $plus.minus$ SD across $N=50$ samples): at $r <= 0.8$, all methods maintain minimal error ($2.5 times 10^(-5)$ to $1.8 times 10^(-4)$); at $r = 0.99$, regularized and Bayesian methods contract collinear siblings toward their centroid, trading small shrinkage bias for topological stability. Dashed slate line confirms broad lineage MSE ($approx 5.4 times 10^(-4)$) is invariant across methods. (C) Intra-lineage sibling cross-talk correlation ($"Corr"(hat(theta)_a, hat(theta)_b)$): BayesPrism and InstaPrism scale smoothly toward $+1.00$, RegDeconv eliminates anti-correlation and jumps to $+0.44$ at $r = 0.99$, whereas NNLS ($-0.130$), Rectangle ($0.000$), CIBERSORT (reimpl., $-0.135$), and CIBERSORTx Docker ($-0.132$) remain trapped in negative cross-talk or zero correlation. (D) Spurious state dropout rate (% true non-zero states collapsed to boundary $0.0$): NNLS collapses $10.0%$, CIBERSORT reimpl. $9.5%$, and CIBERSORTx Docker $7.5%$ of states to false zeros, whereas RegDeconv suppresses dropouts to $0.5%$, and BayesPrism and InstaPrism maintain $0.0%$ dropouts.]
) <fig-collinearity-benchmark>

#figure(
  image("../figures/deconvolution/deconvolution_results_vs_ground_truth.png", width: 95%),
  caption: [Deconvolution ground-truth recovery under severe collinearity ($r = 0.99$). (A) Scatter plot grid of inferred cell fraction ($hat(theta)$) vs. ground truth ($theta^*$) across all six methods with diagonal $y=x$ dashed reference lines, showing $R^2$, RMSE, and Pearson $r$. (B) True vs. predicted cellular composition for representative synthetic samples (Samples 01 to 06), illustrating how unregularized NNLS collapses sibling state mass to zero while RegDeconv, BayesPrism, and InstaPrism preserve proportional representation across synthetic subsets `L1_State_1` to `L4_State_3`.]
) <fig-ground-truth-comparison>

#figure(
  table(
    columns: (1.1fr, 0.9fr, 1.1fr, 1.1fr, 1.1fr, 1.1fr, 1.1fr, 1.1fr, 0.9fr, 0.9fr),
    inset: 4.5pt,
    align: center,
    [*Target $r$*], [*Emp $r$*], [*Unreg $kappa$*], [*Reg $kappa$*], [*NNLS MSE*], [*Reg MSE*], [*NNLS $r_("sib")$*], [*Reg $r_("sib")$*], [*NNLS Drop*], [*Reg Drop*],
    [0.00], [-0.011], [$3.3 times 10^1$], [$3.3 times 10^1$], [0.00002], [0.00002], [-0.122], [-0.121], [0.2%], [0.3%],
    [0.20], [+0.212], [$7.1 times 10^1$], [$7.1 times 10^1$], [0.00004], [0.00004], [-0.159], [-0.161], [0.0%], [0.0%],
    [0.60], [+0.614], [$1.6 times 10^2$], [$1.6 times 10^2$], [0.00008], [0.00009], [-0.094], [-0.090], [1.0%], [0.7%],
    [0.80], [+0.805], [$2.9 times 10^2$], [$2.9 times 10^2$], [0.00017], [0.00018], [-0.081], [-0.103], [0.7%], [0.8%],
    [0.95], [+0.948], [$8.6 times 10^2$], [$8.6 times 10^2$], [0.00033], [0.00051], [-0.094], [-0.087], [4.2%], [2.7%],
    [0.98], [+0.978], [$1.8 times 10^3$], [$1.8 times 10^3$], [0.00077], [0.00073], [+0.004], [+0.021], [7.0%], [5.2%],
    [0.99], [+0.990], [$3.8 times 10^3$], [*$1.6 times 10^3$*], [0.00132], [*0.00290*], [-0.130], [*+0.444*], [10.0%], [*0.5%*]
  ),
  caption: [Continuous synthetic benchmark across the correlation spectrum ($r in [0.0, 0.99]$) comparing unregularized NNLS against Adaptive Deficit RegDeconv ($S=12$ states, $N=50$ samples per correlation level).]
) <tab-collinearity-benchmark>

This continuous-spectrum benchmark establishes three fundamental computational insights:
1. *Zero Overregularization in Well-Conditioned Regimes*: In uncollinear and moderate regimes ($r <= 0.8$, $kappa <= 300$), Adaptive Deficit Regularization injects zero excess curvature ($lambda_("lap") approx 0$), ensuring that state MSE ($2.5 times 10^(-5)$ to $1.8 times 10^(-4)$) is mathematically identical between NNLS and RegDeconv. This resolves the midpoint shrinkage bias observed in fixed-regularization schemes.
2. *Adaptive Curvature Injection at Severe Collinearity*: As sibling state profiles reach extreme collinearity ($r = 0.99$), unregularized Hessian condition numbers surge to $3.8 times 10^3$. The adaptive deficit mechanism bounds $kappa(bold(H)) <= 2,000$ ($kappa = 1,618$), stabilizing the estimation manifold (@fig-collinearity-benchmark a).
3. *Suppression of Boundary Dropouts and Negative Cross-Talk*: Under standard NNLS, ill-conditioning forces $10.0%$ of true non-zero cell states to collapse to numerical boundary zeros ($hat(theta) = 0.0$), with severe negative sibling cross-talk ($r_("sibling") = -0.130$). RegDeconv suppresses false dropouts down to $0.5%$ and restores coherent positive sibling correlation ($r_("sibling") = +0.444$) (@fig-collinearity-benchmark c, d). While smoothing collinear siblings introduces a small centroid shrinkage bias when true proportions are heterogeneous ($"MSE" = 0.0029$ vs $0.0013$ for NNLS), RegDeconv retains less than half the shrinkage error of InstaPrism ($0.0057$) and BayesPrism ($0.0060$) while fully safeguarding manifold topology.

== 3.15 Estimator Bias-Variance Tradeoff and True Biological Dropout Sensitivity
To rigorously dissect why Bayesian and EM hierarchical deconvolution engines (BayesPrism and InstaPrism) exhibit lower fine-state Pearson correlation with ground truth ($r approx 0.53 - 0.55$) at severe collinearity ($r = 0.99$) despite exceptional lineage-level fidelity, we performed two targeted diagnostic stress tests (@fig-variance-dropout, @tab-variance-dropout):

#figure(
  image("../figures/deconvolution/variance_and_dropout_diagnostics.png", width: 95%),
  caption: [Estimator bias-variance decomposition and biological zero dropout stress-testing. (A1) Exact MSE decomposition ($"MSE" = "Bias"^2 + "Variance"$) across $B=30$ independent sequencing draws from the same biological mixture at $r = 0.99$: unregularized NNLS is $95.6%$ estimator variance, whereas BayesPrism and InstaPrism are $>99.5%$ centroid shrinkage bias with near-zero variance. (A2) Estimator variance $"Var"(hat(theta))$ on $log_(10)$ scale: BayesPrism ($2.8 times 10^(-5)$) and InstaPrism ($3.0 times 10^(-6)$) reduce technical noise variance by $45times$ to $400times$ relative to NNLS ($1.2 times 10^(-3)$). (B1) False positive phantom detection mass (%) under true biological zero absence ($theta^* = 0$): Dirichlet shrinkage priors leak $8.3%$ mass into absent sibling states when their lineage is active, and BayesPrism accumulates $15.7%$ total phantom mass across absent lineages. (B2) Reconstruction RMSE on active non-zero subsets.]
) <fig-variance-dropout>

#figure(
  table(
    columns: (1.5fr, 1.1fr, 1.1fr, 1.1fr, 1.1fr, 1.1fr),
    inset: 4.5pt,
    align: center,
    [*Method*], [*Estimator Var*], [*Var % of MSE*], [*Squared Bias*], [*Phantom Sibling*], [*Total Phantom*],
    [Unregularized (NNLS)], [$1.25 times 10^(-3)$], [95.6%], [$9.9 times 10^(-5)$], [1.84%], [2.08%],
    [Rectangle (DWLS-QP)], [$1.02 times 10^(-3)$], [91.4%], [$1.3 times 10^(-4)$], [1.68%], [1.92%],
    [CIBERSORT (reimpl.)], [$1.37 times 10^(-3)$], [89.0%], [$2.1 times 10^(-4)$], [2.26%], [5.69%],
    [RegDeconv (Graph Lap)], [*$5.35 times 10^(-4)$*], [*21.8%*], [$1.94 times 10^(-3)$], [5.64%], [5.84%],
    [InstaPrism], [*$3.0 times 10^(-6)$*], [*0.06%*], [$5.30 times 10^(-3)$], [8.25%], [11.57%],
    [BayesPrism (Gibbs)], [*$2.8 times 10^(-5)$*], [*0.50%*], [$5.65 times 10^(-3)$], [8.35%], [15.74%]
  ),
  caption: [Empirical bias-variance decomposition across $B=30$ sequencing replicates from a fixed tissue sample ($r = 0.99$) and false positive phantom mass on true biological zeros ($N=20$ samples).]
) <tab-variance-dropout>

These diagnostic experiments uncover two key statistical realities:
1. *The Low-Variance Advantage of Hierarchical Priors*: Across repeated technical sequencing draws from the exact same tissue sample, unregularized NNLS, Rectangle, and CIBERSORT fluctuate wildly ($"Var" approx 1.0 - 1.4 times 10^(-3)$), with variance accounting for $89% - 96%$ of total error. In sharp contrast, InstaPrism and BayesPrism achieve *near-zero estimator variance* ($3.0 times 10^(-6)$ and $2.8 times 10^(-5)$), suppressing technical noise variance by up to $400times$. Their lower fine-state Pearson correlation reflects a deterministic centroid contraction ($omega_(t, s) approx 1/3$) rather than instability: they intentionally sacrifice within-lineage fine-scale resolution to eliminate technical replicate variance. RegDeconv provides an intermediate balance, reducing NNLS variance by $57%$ while incurring less than half the shrinkage bias of BayesPrism.
2. *Dirichlet Prior Phantom Leakage on True Zeros*: Under true biological absence ($theta^* = 0$), Dirichlet priors with pseudo-counts necessarily distribute baseline probability mass among all states within an active cell type. Consequently, when an individual cell state is absent, BayesPrism and InstaPrism erroneously allocate $8.3% - 8.4%$ false positive fraction to the missing state, and BayesPrism assigns $7.39%$ to an entirely absent lineage ($15.74%$ total phantom mass). Conversely, sparsity-compatible estimators (NNLS and Rectangle) achieve clean zero preservation ($< 2.1%$ total phantom detection). RegDeconv fully suppresses absent lineages ($0.20%$ phantom mass) while exhibiting modest smoothing leakage ($5.64%$) among active sibling states.

== 3.16 Regimes of Bayesian Superiority: Depth, Reference Drift, Overdispersion, and Phenotypic Plasticity
The preceding bias-variance decomposition demonstrates that hierarchical Bayesian priors (BayesPrism and InstaPrism) and Graph-Laplacian regularization (RegDeconv) introduce deterministic centroid shrinkage to extinguish ill-conditioned estimator variance. This raises a crucial question: under what concrete biological and technical conditions do Bayesian estimators achieve decisively lower total error (MSE) than standard unregularized or quadratic/support-vector deconvolution? We identified and benchmarked four widespread experimental regimes where BayesPrism and InstaPrism systematically outperform NNLS, Rectangle, and CIBERSORT (@fig-superior-regimes, @tab-superior-regimes):

#figure(
  image("../figures/deconvolution/bayesprism_superior_regimes.png", width: 95%),
  caption: [Empirical regimes of Bayesian deconvolution dominance across biological and technical stress tests under severe collinearity ($r = 0.99$). (A) Low sequencing depth ($n_("total") in [1,000, 80,000]$): as read depth drops below $10"k"$, Poisson sampling noise detonates the variance of NNLS ($"MSE" = 0.0152$), Rectangle ($"MSE" = 0.0139$), and CIBERSORT ($"MSE" = 0.0133$), while InstaPrism ($"MSE" = 0.0045$) and BayesPrism ($"MSE" = 0.0066$) achieve up to $3.4times$ lower error. (B) Patient-specific reference expression drift ($sigma_("drift") in [0.0, 0.75]$): multiplicative lognormal distortion of tumor/stromal signatures causes linear least-squares error to surge to $0.0113$, whereas multinomial likelihood dampens gene expression outliers, keeping InstaPrism ($0.0048$) and BayesPrism ($0.0051$) completely flat ($2.4times$ error reduction). (C) Biological overdispersion and transcriptional bursting ($alpha_("disp") in [0.0, 0.30]$): under Negative Binomial extra-Poisson variance, $L_2$-loss penalties over-leverage bursty genes ($"MSE" = 0.0153$), whereas Bayesian count models bound outlier leverage ($"MSE" = 0.0052$, $3.0times$ reduction). (D) Balanced intra-lineage plasticity ($alpha_("within") in [1, 40]$): when collinear sibling states co-occur evenly around parity ($theta_(t,s) approx 1/3$), shrinkage bias drops to zero, and InstaPrism ($"MSE" = 0.00017$, $13.9times$ reduction vs. NNLS) and BayesPrism ($"MSE" = 0.00036$, $6.7times$ reduction) decisively dominate all linear models.]
) <fig-superior-regimes>

#figure(
  table(
    columns: (1.5fr, 1.1fr, 1.1fr, 1.1fr, 1.1fr),
    inset: 4.5pt,
    align: center,
    [*Method*], [*Low Depth ($1"k"$)*], [*Drift ($sigma = 0.75$)*], [*Overdisp. ($alpha = 0.30$)*], [*Plasticity ($alpha = 40$)*],
    [Unregularized (NNLS)], [0.01522], [0.01133], [0.01533], [0.00239],
    [Rectangle (DWLS-QP)], [0.01392], [0.00844], [0.01254], [0.00241],
    [CIBERSORT (reimpl.)], [0.01326], [0.00938], [0.01417], [0.00247],
    [CIBERSORTx (Docker)], [0.01326], [0.00938], [0.01417], [0.00247],
    [RegDeconv (Graph Lap)], [0.01064], [0.00644], [0.01056], [0.00062],
    [BayesPrism (Gibbs)], [*0.00661*], [*0.00513*], [*0.00522*], [*0.00036*],
    [InstaPrism], [*0.00447*], [*0.00480*], [*0.00517*], [*0.00017*]
  ),
  caption: [Benchmark summary of state proportion MSE across four stress regimes at severe collinearity ($r = 0.99$). Bold values highlight top-performing Bayesian estimators.]
) <tab-superior-regimes>

1. *Low Sequencing Depth and Poisson Shot Noise ($n_("total") <= 10,000$)*:
   Under finite bulk sequencing counts, Poisson sampling variance injects severe counting noise into the observed counts $bold(y)$. Because ill-conditioned least squares scales as $"Cov"(hat(bold(theta))) prop 1 / (n_("total") (1 - r))$, unregularized estimators suffer an explosive variance surge at low reads. At $n_("total") = 1,000$ counts, total MSE for unregularized NNLS reaches $0.01522$, Rectangle $0.01392$, and CIBERSORT $0.01326$. In contrast, InstaPrism ($"MSE" = 0.00447$) and BayesPrism ($"MSE" = 0.00661$) remain robust, delivering *over $3.4times$ lower error* (@fig-superior-regimes a, @tab-superior-regimes).

2. *Patient-Specific Reference Drift and Tumor Plasticity ($sigma_("drift") >= 0.50$)*:
   Reference expression profiles $phi_(s,g)$ derived from single-cell atlases inevitably diverge from bulk tumor biopsies due to patient genetics, subclonal evolution, and microenvironmental cross-talk. When individual sample profiles exhibit multiplicative lognormal drift ($phi_(n,s,g) prop phi_(s,g) dot delta_(n,s,g)$ with $log delta ~ cal(N)(0, sigma^2)$), linear $L_2$ regression models suffer severe quadratic penalties from highly drifting genes, driving NNLS MSE up to $0.01133$ and CIBERSORT to $0.00938$ at $sigma = 0.75$. Conversely, the Multinomial log-likelihood $sum_g y_g log sum_s theta_s phi_(s,g)$ naturally normalizes across genes, rendering InstaPrism ($0.00480$) and BayesPrism ($0.00513$) virtually impervious to reference drift (*$2.4times$ lower error*).

3. *Biological Overdispersion and Transcriptional Bursts ($alpha_("disp") >= 0.15$)*:
   Real RNA-seq counts exhibit substantial overdispersion beyond Poisson variance, generated by transcriptional bursting and biological heterogeneity ($y_(n,g) ~ "NegBinomial"(mu_(n,g), alpha)$ with $"Var"(Y) = mu + alpha mu^2$). Outlier burst genes exert outsized leverage on quadratic regression objectives. Under $alpha = 0.30$, NNLS MSE detonates to $0.01533$ and CIBERSORT to $0.01417$. In contrast, InstaPrism ($0.00517$) and BayesPrism ($0.00522$) achieve *over $3.0times$ lower error*, demonstrating that discrete likelihood formulations prevent high-variance outlier genes from distorting proportion estimates.

4. *Balanced Phenotypic Co-occurrence and Plasticity ($alpha_("within") >= 10$)*:
   When collinear sibling states represent continuous phenotypic plasticity within a lineage—such as CD8+ T cell exhaustion states or macrophage activation spectra—they co-exist within the same spatial niche rather than occurring in mutually exclusive, extreme proportions. When sibling state proportions follow a balanced Dirichlet distribution ($alpha_("within") = 40$), the true state proportions co-occur near parity ($theta_(t,s) approx 1/3 theta_t$). Here, the centroid shrinkage prior matches the true data manifold, collapsing shrinkage bias to near-zero while preserving full variance suppression. Under balanced states, *InstaPrism achieves an MSE of $0.000172$* ($13.9times$ lower than NNLS and $14.4times$ lower than CIBERSORT) and *BayesPrism achieves $0.000357$* ($6.7times$ lower than NNLS), while RegDeconv achieves $0.000616$.

== 3.17 Unified 3-Factorial Error Landscape and Algorithm Decision Matrix
To synthesize these findings into a unified, principled taxonomy across the entire computational oncology landscape, we executed a complete 3-factorial in silico benchmark systematically decomposing total error into its exact mathematical constituents ($"MSE" = "Bias"^2 + "Variance"$) across four collinearity levels ($r in [0.0, 0.6, 0.9, 0.99]$), three count parameters ($n_("total") in [2,500, 10,000, 80,000]$), and three biological state architectures (`skewed`, `balanced`, `state_dropout`) across all seven deconvolution tools (@fig-unified-landscape, @fig-dominance-map).

#figure(
  image("../figures/deconvolution/unified_bias_variance_landscape.png", width: 95%),
  caption: [Two-dimensional faceted bias-variance decomposition landscape ($"MSE" = "Bias"^2 + "Variance"$) across collinearity levels $r in [0.0, 0.99]$ (columns) and sequencing count parameters $n_("total") in [2,500, 80,000]$ (rows) across all seven deconvolution methods under heterogeneous biological proportions. Stacked bars represent total state MSE, partitioned into estimator variance $"Var"(hat(theta))$ (purple) and squared bias $"Bias"^2(hat(theta))$ (gold). At low counts ($n_("total") = 2,500$) and severe collinearity ($r = 0.99$), unregularized NNLS, CIBERSORT, and CIBERSORTx Docker explode into pure estimator variance, whereas InstaPrism and BayesPrism maintain near-zero variance, achieving lower total MSE.]
) <fig-unified-landscape>

#figure(
  image("../figures/deconvolution/algorithm_dominance_matrix.png", width: 95%),
  caption: [Empirical algorithm dominance decision map identifying the winning deconvolution method (lowest total MSE) across collinearity $r$ and count parameter $n_("total")$ for each biological state architecture. (Left) Balanced intra-lineage plasticity: InstaPrism dominates almost universally, with RegDeconv optimal at high counts and moderate collinearity. (Center) Heterogeneous proportions: InstaPrism dominates at low counts ($n_("total") <= 10,000$), while NNLS and RegDeconv share the high-count regime. (Right) Structural state dropout ($theta^* = 0$): RegDeconv dominates across moderate-to-high counts by avoiding Dirichlet phantom leakage, while InstaPrism buffers against Poisson shot noise at low counts.]
) <fig-dominance-map>

This unified analysis yields clear, actionable guidelines for algorithm selection in bulk RNA-seq deconvolution:
1. *Low-Pass and Spatial Transcriptomics ($n_("total") <= 10,000$)*: Hierarchical shrinkage engines (*InstaPrism* and *BayesPrism*) are strictly superior across all collinearity levels and biological architectures, eliminating up to $99%$ of the Poisson sampling variance that cripples unregularized least-squares and support-vector estimators.
2. *Phenotypic Plasticity and Balanced Co-occurrence*: When intra-lineage cell states represent continuous gradients (e.g., exhaustion or activation spectra), shrinkage bias vanishes, rendering *InstaPrism* and *RegDeconv* decisively superior ($3.5times$ to $5.0times$ lower error) even at high sequencing depth.
3. *Presence of True Biological Dropouts ($theta^* = 0$)*: When rare cell states or lineages are truly absent, Dirichlet shrinkage priors leak phantom fraction mass into missing subsets ($8% - 16%$). Under these conditions with adequate sequencing depth ($n_("total") >= 10,000$), *RegDeconv* delivers the optimal compromise: bounding ill-conditioned matrix inversion curvature without Dirichlet phantom leakage.

