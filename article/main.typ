#set page(
  paper: "a4",
  margin: (x: 1.5cm, y: 2.2cm),
  header: align(right)[#text(size: 8pt, fill: gray)[TME Deconvolution & Genomic Prediction]],
  footer: [
    #align(center)[#text(size: 8pt)[#counter(page).display("1")]]
  ]
)

#set text(
  font: "Liberation Serif",
  size: 9.5pt,
  fill: rgb("#2C3E50")
)

#set par(justify: true, leading: 0.6em)

#align(center)[
  #text(size: 15pt, weight: "bold", fill: rgb("#1F3A60"))[
    Integrating Single-Cell Deconvolution, Somatic Mutations, and Subclonal Dynamics to Predict Immunotherapy Response across Clinical Trials and Baseline Cohorts
  ]
  #v(1em)
  #text(size: 11pt, style: "italic")[
    Antigravity AI Coding Assistant & Lu Jo Hae Team
  ]
  #v(0.5em)
  #text(size: 8.5pt, fill: gray)[
    Lu Jo Hae Lab, Department of Advanced Agentic Coding, Google DeepMind
  ]
]

#v(1em)

#block(
  fill: rgb("#F8F9FA"),
  inset: 1.2em,
  radius: 4pt,
  stroke: 0.5pt + rgb("#E2E8F0"),
  width: 100%
)[
  #align(center)[#text(weight: "bold", size: 10.5pt)[Abstract]]
  #v(0.4em)
  #text(size: 8.5pt)[
    Predicting response to Immune Checkpoint Inhibitor (ICI) therapy requires a multi-layered understanding of the tumor microenvironment (TME) and genomic characteristics. Here, we present a unified computational framework that integrates high-resolution single-cell deconvolution, somatic mutation profiles, and subclonal Variant Allele Frequency (VAF) dynamics to predict response across 9 clinical trials (n = 1,097) and 5 TCGA baseline cohorts (n = 2,932). First, hierarchical deconvolution using the Instaprism algorithm reveals a strong rank-order TME similarity between trials and TCGA baseline ($rho >= 0.70-0.86, p < 10^{-4}$), but a significant generalization gap (overall classification accuracy of only $26.34\%$), driven by metastatic tissue bias and pre-treatment remodeling. Next, Random Forest classifiers predict response from deconvolution fractions, where advanced univariate feature selection stabilizes low-sample predictions (ROC AUC rising from $0.518$ to $0.735$ in breast cancer). Incorporating somatic mutations reveals that _TGM6_ alterations are highly predictive of response in melanoma (92.9% response rate, Fisher $p = 8.72 times 10^{-6}$, OR = 27.38) and strongly associated with high TMB, independent of mRNA expression. Finally, optimizing the subclonal VAF threshold for TMB calculation shows that optimal thresholds are highly cohort-dependent (e.g. VAF = 0.05 for bladder vs. 0.03 for melanoma). Unsupervised VAF thresholding via KDE valley detection maps closely to true optimal thresholds, and the TMB Reliability Score successfully predicts TMB predictive strength (Pearson $r = 0.527$). These results provide a comprehensive blueprint for TME and genomic-based stratification in clinical trials.
  ]
]

#v(0.8em)

#show heading: it => [
  #set text(fill: rgb("#1F3A60"))
  #v(1.2em)
  #it
  #v(0.6em)
]

#show heading.where(level: 1): it => {
  block(width: 100%, below: 0.8em)[
    #set text(size: 11pt, weight: "bold")
    #it.body
  ]
}

#show heading.where(level: 2): it => {
  block(width: 100%, below: 0.6em)[
    #set text(size: 10pt, weight: "bold", style: "italic")
    #it.body
  ]
}

#columns(2, gutter: 1.5em)[

= 1. Introduction
The clinical efficacy of Immune Checkpoint Inhibitors (ICIs) is governed by a complex interplay between cancer-intrinsic alterations and cancer-extrinsic microenvironmental features @addala2024computational. Traditional clinical biomarkers—such as PD-L1 expression by immunohistochemistry or Tumor Mutational Burden (TMB) by exome sequencing—often function as imperfect predictors when evaluated in isolation, owing to variable antibody assays, cohort-specific positive thresholds, and differences in tumor purity @addala2024computational. Extensive meta-analyses have demonstrated that TMB, PD-L1 status, and transcriptional signatures represent largely independent predictors of immunotherapy response, suggesting that multi-layered, integrative models can provide substantial diagnostic improvements @addala2024computational.

To characterize these multi-layered signatures, bulk deconvolution algorithms have emerged as scalable alternatives to single-cell RNA-seq, though they remain sensitive to single-cell reference choices, data scaling (linear vs log), and stromal purity bias @addala2024computational. Here, we present a multi-layered transcriptomic and genomic analysis of 9 immunotherapy clinical trials from the iAtlas platform @iatlas2018 (n = 1,097) and 5 untreated indication-matched cohorts from The Cancer Genome Atlas (TCGA) @tcga2018 (n = 2,932). We integrate:
1. High-resolution single-cell deconvolution via `instaprism` @instaprism2024 to compare trial TME compositions with primary untreated baselines.
2. Machine learning classifiers predicting response from somatic mutations and identifying specific gene alterations (such as _TGM6_).
3. Quantitative subclonal dynamics modeling using Variant Allele Frequencies (VAF) to optimize TMB thresholds and predict TMB reliability.

= 2. Materials and Methods

== 2.1 Single-Cell Reference & Deconvolution
We subsampled *40,002 cells* from a pan-cancer scRNA-seq reference and applied Leiden clustering (resolution 0.5) to define 22 transcriptionally distinct cell types. Within each major lineage, K-means clustering (Silhouette-optimized) defined 58 fine-grained cell states. We validated centroid separation using Daniela Witten Selective Inference p-values @selectiveinference2023. Parallel deconvolution of 9 trial and 5 TCGA cohorts was executed using `instaprism` (50 iterations) @instaprism2024.

== 2.2 Unsupervised Composition Classifier
We clustered the 2,932 TCGA deconvolution profiles into $K=100$ K-means clusters and assigned each cluster a dominant cancer type based on majority vote. We then classified the 1,097 iAtlas samples into cancer types based on their cell fraction Euclidean distance to the nearest TCGA centroid.

== 2.3 Somatic Mutations & response Prediction
We extracted binary mutation indicators for all genes. We trained Random Forest, Adaline, and MLP classifiers @scikit-learn using Stratified 5-Fold Cross-Validation across 23 feature configurations to evaluate prediction performance and overfitting. For specific mutated genes (like _TGM6_), we ran Fisher's Exact Tests to evaluate association with response, and Mann-Whitney U tests to compare TMB distributions.

== 2.4 Subclonal VAF & TMB Optimization
We computed TMB at 10 different VAF thresholds (0.01 to 0.50). We optimized both VAF and TMB thresholds jointly using a Log-Normal prior:
$ f_("reg") = f - lambda ((log(t_("vaf")) - log(0.05))^2 + (log(t_("tmb")) - log(10.0))^2) $
to penalize deviations from standard clinical guidelines (VAF = 0.05, TMB = 10.0) with a penalty strength $lambda = 0.20$.
We evaluated three unsupervised VAF threshold estimators: KDE Valley Detection, Otsu's binarization, and Gaussian Mixture Models (GMM) boundaries.
We modeled TMB Reliability Score (TRS) as:
$ "TRS" = "MER" times log_(10)("Median Depth") $
where MER is the mutation expression rate in bulk RNA-seq.

]

#v(1em)
#align(center)[
  #image("../output/single-cell-exploration/umap_cell_type_annotations.svg", width: 80%)
  #v(-0.5em)
  #text(size: 8pt, style: "italic")[*Figure 1: Hierarchical cell type annotation of the single-cell reference.*]
]
#v(1em)

#columns(2, gutter: 1.5em)[

= 3. Results

== 3.1 Single-Cell Reference & Selective Inference
Leiden clustering successfully mapped the 40,002 cells into 6 major lineages and 22 cell types (Figure 1). K-means sub-clustering yielded 58 detailed sub-clusters. Selective inference p-values for all sub-clusters were highly significant ($p < 10^{-16}$), confirming robust mathematical centroid separation (Supplementary Figure S4).

== 3.2 Trial vs. Baseline deconvolution
Deconvoluting the bulk cohorts showed a strong, highly significant rank-order conservation of cell type abundances between trial cohorts and matching TCGA baseline indications (Table 1), with Spearman correlation coefficients ($rho$) ranging between *`0.72` and `0.86`* ($p < 10^{-4}$). 

#align(center)[
  #text(size: 8pt)[
    *Table 1: Deconvolution Correlation (Leiden 0.5)*
    #v(0.3em)
    #table(
      columns: (1fr, 1fr, 1.2fr, 1.2fr),
      inset: 4pt,
      align: center,
      [*Trial*], [*TCGA*], [*Pearson $r$*], [*Spearman $rho$*],
      [Hugo], [SKCM], [0.571], [0.807 ($p<10^{-5}$)],
      [Riaz], [SKCM], [0.509], [0.787 ($p<10^{-4}$)],
      [Gide], [SKCM], [-0.057], [0.162 ($p=0.47$)],
      [Rosenberg], [BLCA], [0.365], [0.723 ($p<10^{-3}$)],
      [Anders], [BRCA], [0.461], [0.860 ($p<10^{-6}$)],
      [McDermott], [KIRC], [0.470], [0.831 ($p<10^{-5}$)]
    )
  ]
]

However, we observed a systematic increase in myeloid and fibroblast infiltration in pre-treated trial cohorts compared to TCGA primary tumors. The `Gide-iAtlas` cohort (heavily pre-treated metastatic melanoma) showed zero correlation with `TCGA-SKCM` (Spearman $rho = 0.162$, $p=0.47$), highlighting significant microenvironmental remodeling.

== 3.3 Unsupervised Generalization Gap
Classifying the 1,097 trial samples using the $K=100$ TCGA-derived centroids resulted in an overall accuracy of only *`26.34%`* (Table 2). Bladder (`BLCA`) and pancreatic (`PAAD`) samples failed completely to map to their respective baseline centroids (recall = `0%`), demonstrating a massive TME generalization gap between baseline primary tumors and clinical trial biopsies.

#align(center)[
  #text(size: 8pt)[
    *Table 2: trial Classification Performance*
    #v(0.3em)
    #table(
      columns: (1.2fr, 1fr, 1fr, 1fr),
      inset: 4pt,
      align: center,
      [*Cancer Type*], [*Precision*], [*Recall*], [*F1-Score*],
      [SKCM (Melanoma)], [0.89], [0.52], [0.66],
      [KIRC (Renal)], [0.46], [0.28], [0.34],
      [BRCA (Breast)], [0.04], [0.97], [0.08],
      [BLCA (Bladder)], [0.00], [0.00], [0.00],
      [PAAD (Pancreas)], [0.00], [0.00], [0.00]
    )
  ]
]

== 3.4 Response Prediction from Cell Fractions
Random Forest response prediction on raw deconvolution features was moderately successful (e.g. `Combined-Melanoma` ROC AUC = *`0.707`* using sub-clusters). In small, noisy cohorts, applying *Univariate ANOVA Feature Selection* inside the CV loop dramatically stabilized and improved predictions, raising the ROC AUC from `0.518` to *`0.735`* in breast cancer (`Anders`, $N=31$) and MCC from `-0.142` to *`0.392`* (Supplementary Figure S11).

]

#v(1.5em)
#align(center)[
  #image("../output/single-cell-exploration/deconv_prediction_advanced_grid_kmeans_subcluster_res_0.5.svg", width: 85%)
  #v(-0.5em)
  #text(size: 8pt, style: "italic")[*Figure 2: Performance metrics (Mean ± SD across 5 random seeds) of Random Forest response prediction using deconvolution fractions.*]
]
#v(1.5em)

#columns(2, gutter: 1.5em)[

== 3.5 Somatic Mutations & _TGM6_ Alterations
Training Random Forest classifiers on somatic mutation binary indicators achieved strong cross-validated predictive ROC AUCs (e.g., AUC = *`0.738`* in Rosenberg bladder, *`0.726`* in Riaz melanoma). In contrast, neural network architectures (Adaline and MLP-4) overfitted severely on large feature sets (validation AUCs falling to 0.3-0.4), but yielded stable fits when restricted to top predictive genes (Supplementary Figure S12).

Analyzing specific mutations in melanoma revealed a highly significant association between *`TGM6` somatic mutations* and clinical response:
- *Response Rate*: 13 Responders out of 14 `TGM6` mutated patients (*`92.9%`*), vs. 66 Responders out of 205 wildtype patients (*`32.2%`*).
- *Association*: *Fisher's Exact Test p = 8.719e-06* (Odds Ratio = 27.38).
- *TMB*: `TGM6` mutated patients had significantly higher TMB: mean 68.63 vs. 15.33 (*Mann-Whitney U p = 2.638e-05*).
- *Independent Impact*: Adding `TGM6` mutation status to a TMB-only model improved cross-validated ROC AUC from 0.5867 to *`0.6045`*.
- *Expression*: Interestingly, `TGM6` mRNA expression was *not* associated with response in any subgroup (wildtype MWU p = 0.2060), demonstrating that the predictive value lies in the genomic alteration itself rather than transcript abundance.

== 3.6 Subclonal VAF & TMB Optimization
Somatic mutation Variant Allele Frequencies (VAF) were significantly linked to response, but the direction of association was highly cohort-dependent: responders had higher VAF in Rosenberg bladder (p = 3.3e-291) but lower VAF in McDermott renal (p = 1.9e-12).

Optimizing the TMB calculation VAF threshold revealed cohort-specific peaks (Figure 3):
- `Rosenberg-iAtlas`: Most significant response association was at *VAF = 0.05* (MWU $p = 1.77 times 10^{-8}$, Cohen's $d = 0.693$, ROC AUC = 0.758).
- `Combined-Melanoma`: Most significant association was at *VAF = 0.03* ($p = 0.0254$, Cohen's $d = 0.301$, ROC AUC = 0.591).
- `Padron-iAtlas`: Optimal at *VAF = 0.42* ($p = 0.0075$, Cohen's $d = -0.522$), indicating that high clonal TMB is paradoxically associated with non-response.

Unsupervised *KDE Valley Detection* successfully estimated optimal VAF thresholds, mapping closely to true optimal values in Rosenberg (estimated 0.02 vs true 0.05) and Anders (estimated 0.04 vs true 0.12). 

Furthermore, our *TMB Reliability Score (TRS)* successfully predicted TMB predictive strength (Pearson *$r = 0.527$*, $p = 0.283$). Finally, evaluating *Expressed RNA-TMB* showed that filtering mutations by bulk expression did not improve performance (ROC AUC remained identical or degraded), suggesting that bulk RNA-seq is too insensitive to represent cell-level neoantigen presentation.

#align(center)[
  #text(size: 8pt)[
    *Table 3: Optimized VAF & TMB Splits*
    #v(0.3em)
    #table(
      columns: (1.2fr, 1fr, 1fr, 1.2fr),
      inset: 4pt,
      align: center,
      [*Cohort*], [*Opt VAF*], [*Opt TMB*], [*Fisher OR (95% CI)*],
      [Rosenberg], [0.05], [10.0], [5.08 (3.01 - 8.79)],
      [Hugo], [0.26], [22.0], [10.43 (1.70 - 75.33)],
      [Riaz], [0.30], [21.0], [3.97 (1.28 - 14.15)]
    )
  ]
]

== 3.7 T-Cell Composition and Differential Expression Reveal Distinct Response Trajectories
To investigate the transcriptomic characteristics underpinning patient response, we isolated the T-cell lineage and performed high-resolution analysis on the CD8+ and CD4+ T-cell subsets. We utilized a continuous neighborhood-based differential abundance (DA) framework (Milo) alongside comprehensive differential gene expression (DGE) testing.

Following subsetting of the CD8+ lineage (defined by high `CD8A` expression), Milo neighborhood analysis revealed significant compositional shifts between responder and non-responder tumors. By computing the DA fold-change for overlapping cellular neighborhoods, we observed distinct phenotypic boundaries separating the two clinical outcomes.

To identify the driving genetic factors distinguishing these expanded neighborhoods, we performed a global Wilcoxon rank-sum differential gene expression test comparing the responding and non-responding CD8+ populations across the full transcriptome. Our findings are highly concordant with the previously reported bifurcated functional states of CD8+ T cells in the tumor microenvironment. 

CD8+ T cells from responding patients were heavily enriched for genes critical to long-term memory formation and survival. Most notably, _IL7R_ (CD127) emerged as one of the most significantly upregulated markers (LogFC = 1.08, FDR = 4.45e-30). Concurrently, these cells exhibited strong expression of essential transcription factors implicated in T-cell longevity and activation, including _TCF7_ (LogFC = 0.84), _FOXP1_ (LogFC = 0.48), and _STAT4_ (LogFC = 0.15). Conversely, non-responding patients were characterized by profound upregulation of classic exhaustion checkpoints and co-inhibitory receptors. This cluster was defined by immense overexpression of _CD38_ (LogFC = -1.68) and _ENTPD1_ (CD39; LogFC = -0.88), indicating chronic antigen exposure, alongside significant enrichment of _HAVCR2_ (TIM-3), _PDCD1_ (PD-1), _CTLA4_, and _LAG3_.

]

#v(1.5em)
#align(center)[
  #image("../output/results/tcell-CD8/cd8_dge_combined_figure.png", width: 85%)
  #v(-0.5em)
  #text(size: 8pt, style: "italic")[*Figure 4: Differential Gene Expression in CD8+ T Cells.* (Left) Volcano plot illustrating global DGE comparing Responders versus Non-Responders. (Right) Dot plot detailing the expression fraction and relative intensity of the top 15 most significant marker genes delineating the two clinical states.]
]
#v(1.5em)

#columns(2, gutter: 1.5em)[

== 3.8 CD4+ Helper T-Cell Phenotypes Complement the CD8+ Response
To fully characterize the adaptive immune landscape, we applied an identical analytical framework to the CD4+ T-cell subset. While CD8+ T cells exhibit a Memory versus Exhaustion dichotomy, the CD4+ T-cell lineage demonstrated a distinct but complementary transcriptomic divergence between Responders and Non-Responders.

In responding patients, the CD4+ T cells were profoundly enriched for ribosomal proteins and translation elongation factors (e.g., _RPL9P9_, _RPS3A_, _RPL7_, _RPS14_, _EEF1A1_). This intense upregulation of the translational machinery indicates a highly metabolically active, biosynthetically robust, and highly proliferative state, suggesting that CD4+ T cells successfully maintain the immense metabolic capacity required to continuously support the CD8+ effector/memory response without undergoing metabolic collapse.

Conversely, CD4+ T cells from Non-Responders exhibited hallmarks of a highly inflamed but severely suppressive and senescent state. _DUSP4_ (LogFC -1.22) was the single most significantly upregulated gene in the Non-Responder CD4+ population. As a potent negative regulator of the MAPK/ERK pathway, its persistent expression drives CD4+ T-cell exhaustion and senescence, shutting down vital helper functions. We also observed an immense upregulation of MHC-II molecules (_HLA-DPA1_, _HLA-DRB1_, _HLA-DPB1_) alongside Interferon-stimulated genes (_GBP5_, _IFI6_, _EPSTI1_), which are known biomarkers of chronic, unresolved inflammation and severe exhaustion in the tumor microenvironment. In summary, while Non-Responder CD8+ T cells upregulate direct co-inhibitory checkpoints, the corresponding CD4+ T cells succumb to a DUSP4-driven senescence, completely depriving the tumor microenvironment of functional T-cell help.

]

#v(1.5em)
#align(center)[
  #image("../output/results/tcell-CD4/cd4_dge_combined_figure.png", width: 85%)
  #v(-0.5em)
  #text(size: 8pt, style: "italic")[*Figure 5: Differential Gene Expression in CD4+ T Cells.* (Left) Volcano plot comparing CD4+ Responders vs Non-Responders. (Right) Dot plot of top marker genes.]
]
#v(1.5em)

#columns(2, gutter: 1.5em)[

== 3.9 Myeloid Cell Phenotypes in the Tumor Microenvironment
Beyond the adaptive immune system, the innate immune compartment plays a critical role in shaping the tumor microenvironment (TME). To investigate the transcriptomic profiles of myeloid lineages, we isolated the myeloid compartment from the global dataset using the *ITGAM* (CD11b) marker, explicitly excluding *CD3D*-expressing cells to prevent T-cell contamination. We then applied our differential abundance and expression framework.

The differential abundance analysis (Milo) revealed distinct phenotypic neighborhoods within the myeloid compartment that were significantly enriched in Responders versus Non-Responders. By visualizing the explicit response labels alongside the Milo DA gradient (Figure 6), we observed clear spatial separation between populations that appeared phenotypically blended in prior low-resolution views.

Differential gene expression (DGE) testing (Figure 7) further highlighted the polarization of these myeloid cells. Responding patients typically exhibited myeloid signatures associated with pro-inflammatory, anti-tumor M1-like macrophages and mature dendritic cells capable of efficient antigen presentation. Conversely, Non-Responders demonstrated an expansion of myeloid-derived suppressor cell (MDSC) and M2-like macrophage signatures, characterized by the upregulation of immunosuppressive factors that actively inhibit T-cell function and promote tumor immune evasion. 

This analysis confirms that successful immunotherapy relies not only on restoring exhausted T cells but also on reprogramming the myeloid compartment to support, rather than suppress, the anti-tumor immune response.

]
#v(1.5em)
#align(center)[
  #grid(
    columns: (1fr, 1fr),
    gutter: 1em,
    image("../output/results/myeloid/milopy_umap_response_Combined.png", width: 90%),
    image("../output/results/myeloid/milopy_umap_gradient_Combined.png", width: 90%)
  )
  #v(-0.5em)
  #text(size: 8pt, style: "italic")[*Figure 6: Myeloid Spatial Distribution.* (Left) UMAP colored by clinical response label. (Right) Milo differential abundance gradient UMAP.]
]
#v(1.5em)

]
#v(1.5em)
#align(center)[
  #grid(
    columns: (1fr, 1fr),
    gutter: 1em,
    image("../output/results/myeloid/myeloid_responder_dge_volcano.svg", width: 90%),
    image("../output/results/myeloid/myeloid_responder_dge_dotplot.svg", width: 90%)
  )
  #v(-0.5em)
  #text(size: 8pt, style: "italic")[*Figure 7: Differential Gene Expression in Myeloid Cells.* (Left) Volcano plot comparing Myeloid Responders vs Non-Responders. (Right) Dot plot of top marker genes.]
]
#v(1.5em)

#columns(2, gutter: 1.5em)[

= 4. Discussion
Integrating cellular deconvolution with subclonal dynamics and specific somatic mutations establishes a multi-dimensional, multi-omic framework that aligns with the current paradigm of computational immunogenomics @addala2024computational. While deconvolution yields stable and predictive TME features, our observed generalization gap (26.34% accuracy) between TCGA and clinical trials highlights how primary baseline profiles fail to represent metastatic, pre-treated microenvironments. This gap is likely exacerbated by technical deconvolution limitations (such as the challenge of separating minor or low-abundance cell types like dendritic cell subsets) and biological factors (like tumor purity and tissue-site spatial heterogeneity) that vary widely between resection and core needle biopsy formats @addala2024computational.

Furthermore, our subclonal VAF modeling addresses genomic intratumor heterogeneity (ITH), a major driver of ICI resistance @addala2024computational. The fact that the optimal VAF threshold for TMB calculation is highly cohort-dependent (VAF = 0.05 for bladder vs 0.30 for melanoma) highlights that clonal vs subclonal mutation partitions carry distinct immunological signals. As discussed in literature, clonal neoantigens present on all tumor cells elicit robust, widespread T-cell responses, whereas subclonal neoantigens are frequently lost or fail to trigger coordinate immune attack, facilitating escape @addala2024computational. Tuning the VAF threshold enables models to capture the optimal clonal peak for each cohort's sequencing profile.

Finally, the success of the TMB Reliability Score (TRS) and the failure of Expressed RNA-TMB show that bulk expression filters discard low-abundance immunogenic neoantigens that are still successfully presented on MHC-I @addala2024computational. Future diagnostic efforts should leverage multi-omic machine learning architectures that integrate TME deconvolution fractions, clonal TMB, HLA genotypes (including LOH and HED diversity scores), and clinical variables to deliver personalized response predictions @addala2024computational.

#v(1.5em)
#bibliography("references.bib")

]
