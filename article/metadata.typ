// Metadata for the article: Title, authors, abstract, keywords

#let paper_title = [
  Integrating Single-Cell Deconvolution, Somatic Mutations, and Subclonal Dynamics to Predict Immunotherapy Response across Clinical Trials and Baseline Cohorts
]

#let paper_authors = [
  Antigravity AI Coding Assistant & Lu Jo Hae Team
]

#let paper_affiliations = [
  Lu Jo Hae Lab, Department of Advanced Agentic Coding, Google DeepMind
]

#let paper_abstract = [
  Predicting response to Immune Checkpoint Inhibitor (ICI) therapy requires a multi-layered understanding of the tumor microenvironment (TME) and genomic characteristics. Here, we present a unified computational framework that integrates single-cell neighborhood differential abundance (Milo), high-resolution pseudobulk deconvolution, somatic mutations, and subclonal Variant Allele Frequency (VAF) dynamics across 9 clinical trials ($n = 1,097$), 5 TCGA baseline cohorts ($n = 2,932$), and matched clinical melanoma biopsies ($n = 51$). We first benchmark bulk deconvolution directly against single-cell graph DA, demonstrating strong global concordance in metastatic melanoma (Spearman $rho = 0.853, p = 4.18 times 10^(-4)$) and validating naive B cells as concordant responder biomarkers ($hat(beta) = +2.066, p = 3.6 times 10^(-4)$) alongside immunosuppressive macrophages as non-responder biomarkers ($hat(beta) = -1.562, p = 0.0039$), while identifying localized discordance in collinear cytotoxic T-cell sub-lineages. Multi-replicate perturbation stress-testing ($N=5$ replicates, 485 simulations) across 8 distortion modes reveals that patient-specific marker dysregulation induces catastrophic deconvolution failure ($rho$ collapsing to $-0.043 plus.minus 0.473$, sign agreement $47.5\%$), while cell size asymmetry drives severe rank collapse ($rho$ falling to $0.252 plus.minus 0.151$). Compound evaluation identifies the Clinical Core Needle Biopsy regime as the most vulnerable setting ($rho = 0.271$, sign agreement $52.5\%$), and a 2D interaction surface confirms that cell size variation drives rank loss while state activation drives directional sign inversion. Combined with predictive somatic mutations (_TGM6_ alterations yielding OR = 27.38, $p = 8.72 times 10^(-6)$) and unsupervised subclonal VAF optimization, these results establish rigorous mathematical and empirical rules for TME biomarker discovery.
]

#let paper_keywords = (
  "Immunotherapy",
  "Tumor Microenvironment",
  "Single-Cell Deconvolution",
  "Differential Abundance",
  "Variant Allele Frequency",
  "Tumor Mutational Burden",
  "Biomarker Discovery"
)
