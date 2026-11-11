// Supplementary captions for per-atlas Consensus-Meta42 comparisons against pan-cancer T-cell and CAF reference atlases
// plus the unified 42-state Triad specificity matrix.
// Synchronized with outputs/figures in output/results/pan_cancer_atlas_comparison/

#let _atlas_caption_tail(compartment, adj_genes, n_rows) = [
  *(b)* Marker-Anchored Hybrid proximity matrix between the #n_rows `Consensus-Meta42` #compartment states (rows) and every substate of the atlas (columns). A strict one-to-one assignment, obtained by solving a linear sum assignment on the hybrid score, is marked by `●`; no Meta42 state and no atlas substate is used more than once, and matched columns are listed first.
  *(c)* `Consensus-Meta42` `Pan-iAtlas-Meta` bulk deconvolution log-odds effects $hat(beta)_("meta")$ ($N = 1,015$). Panel (a) bars give the meta-analytic logistic log-odds $hat(beta)_("ICB")^("adj")$ of each atlas substate signature across $11$ bulk checkpoint blockade cohorts ($N = 1,058$), adjusted within each trial for lineage load (#adj_genes); asterisks denote Benjamini–Hochberg FDR (`*` $q < 0.05$, `**` $q < 0.01$, `***` $q < 10^(-3)$, `****` $q < 10^(-4)$).
]

#let caption_supp_atlas_zheng_tcell = [
  *`Consensus-Meta42` T/NK states aligned one-to-one with the Zheng et al. pan-cancer T-cell atlas.*
  *(a)* All $41$ meta-clusters of `Zheng_TCell_Atlas` ($397,810$ T cells, $70$ cohorts, $21$ cancer types) @zhengPancancerSinglecellLandscape2021 with their pan-T-adjusted ICB association. Strongest assignments: `M06` $<->$ `CD8.Temra` ($+0.560$), `M40` $<->$ `CD4.TNFRSF9+ Treg` ($+0.540$), `M24` $<->$ `CD4.IL21+ Tfh` ($+0.492$), `M04` $<->$ `CD8.ZNF683+CXCR6+ Trm` ($+0.405$).
  #_atlas_caption_tail([T/NK], [`CD3D`, `CD3E`, `CD3G`, `TRAC`, `TRBC1/2`, `CD2`, `LCK`], 16)
]

#let caption_supp_atlas_chu_tcell = [
  *`Consensus-Meta42` T/NK states aligned one-to-one with the Chu et al. pan-cancer T-cell taxonomy.*
  *(a)* All $39$ fine T-cell states of `Chu_TCell_Atlas` @chuPancancerCellAtlas2023, including heat-shock `CD8_c4_Tstr` and `CD4_c4_Tstr`, represented here by curated marker-signature centroids (`ProfileOrigin.CURATED_SIGNATURE_PRIOR`) because the portal Seurat object was not staged locally. Strongest assignments: `M40` $<->$ `CD4_eTreg` ($+0.334$), `M20` $<->$ `CD4_Tcm` ($+0.258$), `M06` $<->$ `CD8_Teff` ($+0.245$), `M25` $<->$ `CD8_Tex` ($+0.241$).
  #_atlas_caption_tail([T/NK], [`CD3D`, `CD3E`, `CD3G`, `TRAC`, `TRBC1/2`, `CD2`, `LCK`], 16)
]

#let caption_supp_atlas_shi_tnk = [
  *`Consensus-Meta42` T/NK states aligned one-to-one with the T/NK compartment of the Shi et al. cross-tissue atlas.*
  *(a)* All $28$ T/NK clusters of `Shi_PanCancer_Atlas` (stratified $138,452$ cells streamed from the $4,146,975$-cell atlas) @shiCrosstissueMulticellularCoordination2025. Strongest assignments: `M06` $<->$ `I01_CD16hiNK_CREM` ($+0.489$), `M40` $<->$ `CD4T10_Th_ISG15` ($+0.489$), `M04` $<->$ `CD8T02_Tem_GZMK` ($+0.381$), `M01` $<->$ `CD4T03_Tn_NR4A1` ($+0.363$).
  #_atlas_caption_tail([T/NK], [`CD3D`, `CD3E`, `CD3G`, `TRAC`, `TRBC1/2`, `CD2`, `LCK`], 16)
]

#let caption_supp_atlas_gao_caf = [
  *`Consensus-Meta42` stromal states aligned one-to-one with the Gao et al. pan-cancer fibroblast atlas.*
  *(a)* All $20$ fibroblast subtypes `c01`–`c20` of `Gao_CAF_Atlas` ($269,899$ cells) @gaoCrosstissueHumanFibroblast2024. Assignments: `M17` $<->$ `c16_SFRP2_SFRP4_matFib` ($+0.292$), `M02` $<->$ `c17_STMN1_MDK_prolifCAF` ($+0.270$), `M26` $<->$ `c12_C3_KRT18_mesoFib` ($+0.215$), `M12` $<->$ `c07_MMP3_CXCL8_iCAF` ($+0.190$), `M16` $<->$ `c08_RGS5_NOTCH3_Pericyte` ($+0.137$). `c20_CXCL12_CCL19_immunoCAF` is the most responder-enriched substate ($q < 0.01$).
  #_atlas_caption_tail([CAF and mesenchymal-context], [`COL1A1`, `COL1A2`, `COL3A1`, `DCN`, `LUM`, `PDGFRA`], 6)
]

#let caption_supp_atlas_cords_caf = [
  *`Consensus-Meta42` stromal states aligned one-to-one with the Cords et al. CAF phenotypes.*
  *(a)* All $11$ tumor CAF and stromal phenotypes of `Cords_CAF_Atlas` ($29,741$ cells) @cordsCancerassociatedFibroblastPhenotypes2024 @cordsCancerassociatedFibroblastClassification2023. Assignments: `M17` $<->$ `mCAF_MMP11_POSTN` ($+0.291$), `M02` $<->$ `ifnCAF_IDO1_CXCL9` ($+0.271$), `M12` $<->$ `hsp_tpCAF` ($+0.227$), `M16` $<->$ `Pericyte_RGS5_PDGFRB` ($+0.169$), `M26` $<->$ `dCAF_MKI67_TOP2A` ($+0.111$). The interferon-licensed `ifnCAF_IDO1_CXCL9` is strongly responder-enriched ($q = 2.2 times 10^(-6)$).
  #_atlas_caption_tail([CAF and mesenchymal-context], [`COL1A1`, `COL1A2`, `COL3A1`, `DCN`, `LUM`, `PDGFRA`], 6)
]

#let caption_supp_atlas_buechler_caf = [
  *`Consensus-Meta42` stromal states aligned one-to-one with the Buechler et al. cross-tissue fibroblast archetypes.*
  *(a)* All $6$ archetypes of `Buechler_CAF_Atlas` ($10,355$ human fibroblasts) @buechlerCrosstissueOrganizationFibroblast2021. Assignments: `M17` $<->$ `0_COL3A1_panCAF` ($+0.204$), `M02` $<->$ `5_CCL19_FRC_TLS` ($+0.094$), `M16` $<->$ `4_NPNT_alveolar` ($+0.079$), `M26` $<->$ `2_PI16_advProgenitor` ($+0.079$), `M23` $<->$ `1_LRRC15_myCAF` ($+0.048$). With only $6$ archetypes, similarities are low in absolute terms.
  #_atlas_caption_tail([CAF and mesenchymal-context], [`COL1A1`, `COL1A2`, `COL3A1`, `DCN`, `LUM`, `PDGFRA`], 6)
]

#let caption_supp_atlas_shi_caf = [
  *`Consensus-Meta42` stromal states aligned one-to-one with the stromal compartment of the Shi et al. cross-tissue atlas.*
  *(a)* All $19$ stromal clusters of `Shi_PanCancer_Atlas` (stratified $70,593$ cells) @shiCrosstissueMulticellularCoordination2025. Assignments: `M02` $<->$ `S12_Fb_MMP11` ($+0.491$), `M17` $<->$ `S04_Fb_CD9` ($+0.439$), `M12` $<->$ `S08_Fb_PTGES` ($+0.303$), `M16` $<->$ `S03_Fb_LAMC1` ($+0.281$), `M26` $<->$ `S07_Fb_TCF21` ($+0.147$). The `S12_Fb_MMP11` signature is dominated by lymphocyte transcripts (`GNLY`, `CCL5`, `NKG7`) and should be read as a likely doublet or ambient-contamination state rather than a fibroblast program.
  #_atlas_caption_tail([CAF and mesenchymal-context], [`COL1A1`, `COL1A2`, `COL3A1`, `DCN`, `LUM`, `PDGFRA`], 6)
]

#let caption_supp_triad_all42_specificity = [
  *Unified $42$-state block-diagonal transcriptomic specificity of `Consensus-Meta42` across external pan-cancer T/NK, Stromal/CAF, and Myeloid/TAM single-cell reference atlases.*
  *(a)* Pan-lineage hybrid affinity matrix evaluating all $42$ `Consensus-Meta42` deconvolution meta-states (ordered by cellular lineage compartment) against empirical single-cell tumor-infiltrating substate centroids spanning the three-compartment TME Triad (`Zheng_TCell_Atlas` @zhengPancancerSinglecellLandscape2021, `Shi_PanCancer_Atlas` @shiCrosstissueMulticellularCoordination2025, `Gao_CAF_Atlas` @gaoCrosstissueHumanFibroblast2024, `Cords_CAF_Atlas` @cordsCancerassociatedFibroblastPhenotypes2024, `Buechler_CAF_Atlas` @buechlerCrosstissueOrganizationFibroblast2021, and `Coulton_TAM_Atlas` @coultonUsingPancancerAtlas2024, centered relative to the whole-tumor $42$-lineage baseline). Of the $31$ states belonging to T/NK, Stromal/CAF, and Myeloid/TAM lineages, $74.2%$ ($23/31$; mean cognate-versus-off-target margin $+0.134$) assign top-1 to their cognate external compartment.
  *(b)* Row-wise peak Triad affinity colored by `Consensus-Meta42` lineage ontology.
]
