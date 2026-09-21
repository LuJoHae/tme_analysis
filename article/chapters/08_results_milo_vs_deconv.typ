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
- *Collinear Discordance*: Strikingly, `01_Tem/Trm cytotoxic T cells` displayed marked discordance in pre-treatment biopsies (Milo $overline("logFC") = +1.726$, deconvolution $hat(beta) = -0.509$). Because 5 of the 12 clusters represented closely related cytotoxic CD8+ subsets, deconvolution subtracted fraction from `01_Tem/Trm` to mathematically offset the depletion of collinear `07_Tem/Trm`.

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
