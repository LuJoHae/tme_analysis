== 3.7 CD8+ T-Cell Lineage: Memory Longevity vs. Chronic Exhaustion
To dissect the functional adaptive states underpinning patient response, we isolated CD8+ T cells (defined by high _CD8A_ expression) and performed high-resolution continuous neighborhood differential abundance (DA) using Milo, coupled with transcriptome-wide Wilcoxon rank-sum differential gene expression (DGE) testing @vanderleunCD8CellStates2020.

Neighborhood DA analysis uncovered marked phenotypic segregation between responders and non-responders. CD8+ T cells from responding tumors were significantly enriched for genes promoting long-term memory formation and cytokine-independent survival (@fig-cd8-dge). In particular, _IL7R_ (CD127) emerged as the single most significantly upregulated marker ($log_2 "FC" = 1.08$, $"FDR" = 4.45 times 10^(-30)$), co-expressed with canonical memory/stem-like transcription factors including _TCF7_ ($log_2 "FC" = 0.84$), _FOXP1_ ($log_2 "FC" = 0.48$), and _STAT4_ ($log_2 "FC" = 0.15$).

#figure(
  image("../figures/immune-subsets/cd8_dge_combined_figure.png", width: 90%),
  caption: [Differential Gene Expression in CD8+ T Cells. (Left) Volcano plot comparing responding versus non-responding patients. (Right) Dot plot of top 15 marker genes characterizing memory persistence versus terminal exhaustion.]
) <fig-cd8-dge>

Conversely, non-responding patients harbored CD8+ T cells characterized by terminal exhaustion and chronic antigen exposure. This state was defined by substantial overexpression of _CD38_ ($log_2 "FC" = -1.68$) and the ectonucleotidase _ENTPD1_ (CD39; $log_2 "FC" = -0.88$), accompanied by co-inhibitory checkpoint receptors including _HAVCR2_ (TIM-3), _PDCD1_ (PD-1), _CTLA4_, and _LAG3_.

== 3.8 CD4+ Helper T-Cell Phenotypes Complement Adaptive Response
Applying the same analytical framework to CD4+ helper T cells revealed a distinct, complementary transcriptomic divergence. While CD8+ T cells segregated along a memory versus exhaustion axis, CD4+ T cells diverged between high biosynthetic proliferation and suppressive senescence (@fig-cd4-dge).

#figure(
  image("../figures/immune-subsets/cd4_dge_combined_figure.png", width: 90%),
  caption: [Differential Gene Expression in CD4+ T Cells. (Left) Volcano plot comparing CD4+ responders versus non-responders. (Right) Dot plot of marker genes illustrating translational machinery vs. DUSP4-mediated senescence.]
) <fig-cd4-dge>

In responders, CD4+ cells displayed strong enrichment for ribosomal proteins and translation elongation factors (_RPL9P9_, _RPS3A_, _RPL7_, _RPS14_, _EEF1A1_), reflecting a highly active metabolic and translational machinery required to sustain durable helper support. In non-responders, CD4+ T cells were dominated by _DUSP4_ upregulation ($log_2 "FC" = -1.22$). As a negative feedback regulator of MAPK/ERK signaling, sustained _DUSP4_ expression induces premature senescence and helper dysfunction. Non-responder CD4+ cells also upregulated MHC-II alleles (_HLA-DPA1_, _HLA-DRB1_) and interferon-stimulated genes (_GBP5_, _IFI6_, _EPSTI1_), hallmarks of unresolved chronic inflammatory stress.

== 3.9 Myeloid Reprogramming in the Tumor Microenvironment
Beyond lymphoid lineages, we analyzed the innate compartment by isolating myeloid cells (_ITGAM_/CD11b positive, _CD3D_ negative). Milo differential abundance analysis resolved distinct cellular neighborhoods significantly polarized between clinical outcomes (@fig-myeloid-spatial).

#figure(
  grid(
    columns: (1fr, 1fr),
    gutter: 1.2em,
    image("../figures/immune-subsets/milopy_umap_response_Combined.png", width: 100%),
    image("../figures/immune-subsets/milopy_umap_gradient_Combined.png", width: 100%)
  ),
  caption: [Myeloid spatial landscape. (Left) UMAP colored by clinical response status. (Right) Milo differential abundance log-fold change gradient resolving pro- vs anti-tumor myeloid neighborhoods.]
) <fig-myeloid-spatial>

Differential expression analysis confirmed that responding tumors were enriched for pro-inflammatory M1-like macrophages and mature dendritic cell signatures capable of antigen presentation (@fig-myeloid-dge). Non-responders showed pervasive accumulation of myeloid-derived suppressor cell (MDSC) and M2-polarized immunosuppressive programs that inhibit T-cell infiltration and survival, demonstrating that successful checkpoint blockade requires coordinated adaptive activation and innate myeloid reprogramming.

#figure(
  grid(
    columns: (1fr, 1fr),
    gutter: 1.2em,
    image("../figures/immune-subsets/myeloid_responder_dge_volcano.png", width: 100%),
    image("../figures/immune-subsets/myeloid_responder_dge_dotplot.png", width: 100%)
  ),
  caption: [Differential gene expression in the myeloid compartment. (Left) Volcano plot displaying M1/M2 polarization. (Right) Dot plot of dominant innate marker genes distinguishing responders from non-responders.]
) <fig-myeloid-dge>
