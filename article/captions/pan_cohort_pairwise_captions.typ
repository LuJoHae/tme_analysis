// Caption module for the pan-cohort Harmony-integrated pairwise (cross-compartment) Milo interaction atlas

#let caption_pan_cohort_core_hubs_main = [
  *Core pan-cancer host TME hub interactions ($13 times 13$ cell states) on the unified nine-cohort `Harmony I2` atlas.*
  Signed $Z_("ACAT")$ (omnibus Cauchy combination of the discrete $2 times 2$ and continuous `log_ratio_balance` / `co_abundance_product` 3-leg `Berger IUT` tests) for the `13` key cell states retaining $>= 2$ verified partners at $"FDR"_("ACAT") < 0.10$, or $|Z_("ACAT")| >= 3.25$ with stratified permutation $p <= 0.01$, clustered within compartments and grouped into immune and stromal/endothelial super-groups. Black borders mark verified pairs ($"FDR"_("ACAT") < 0.10$); stars denote stratified within-`cohort::substudy` permutation significance ($*$ $p < 0.05$, $**$ $p < 0.01$, $***$ $p < 0.005$). The panel captures `19` of `22` (`86.4%`) verified host pairs.
]

#let caption_pan_cohort_unified_umap_main = [
  *Unified nine-cohort `Harmony I2` UMAP architecture of verified cross-compartment pairwise neighborhoods ($N = 613$ biopsies, `20` `CMH` strata).*
  *(a)* Dual-tier Bezier chords between Meta42 state medoids (bold: $"FDR"_("ACAT") < 0.10$; thin: stratified permutation $p < 0.05$ only), coloured and dashed by pairwise regime, with endpoint polarity glyphs.
  *(b)* Emergent *de novo* pairwise hubs that are non-significant in univariate Milo ($"FDR" >= 0.10$).
  *(c)* Seed-conditioned partner $Z_("ACAT")$ field around the lead anchor.
  *(d)* Lead-pair $k$NN spotlight with $2 times 2$ biopsy response inset (Wilson $95%$ CI).
]

#let caption_pan_cohort_pair_grid_main = [
  *Broad cell-type stratified UMAP atlas of all significant cross-compartment interactions.*
  Panels *(a–k)* show the `11` of `15` unordered pairs of the six broad host cell types (`T-NK`, `B`, `Plasma`, `Myeloid/Macrophages`, `CAF`, `Endothelial`) that harbor at least one significant interaction, encompassing `48` chords in total (`22` bold at $"FDR"_("ACAT") < 0.10$ and `26` thin at stratified permutation $p < 0.05$ only; the `4` empty pairs are omitted). Vivid points mark the lineage cells of interacting states, pale points the remaining lineage context. Panel *(l)* gives the chord, regime and polarity key.
]

#let caption_supp_pan_cohort_harmony_audit = [
  *Outcome-blind `harmonypy` integration benchmark across all nine clinical ICB single-cell cohorts (`I0`–`I3`).*
  Cohort, sub-study and donor `iLISI`, compartment and `Meta42` `cLISI`, `kBET` acceptance, cohort entropy and mixing-QC pass rate for unintegrated PCA (`I0`), `Harmony` on cohort (`I1`), cohort + sub-study (`I2`, selected as Pareto-optimal) and cohort + sub-study + donor (`I3`).
]

#let caption_supp_pan_cohort_full_heatmap = [
  *Full $26 times 26$ pan-cancer host TME cross-compartment $Z_("ACAT")$ heatmap.*
  Signed $Z_("ACAT")$ across `213` of `215` evaluated cross-compartment state pairs, hierarchically clustered within compartment and grouped into immune and stromal/endothelial super-groups. Black borders: $"FDR"_("ACAT") < 0.10$ (`22` pairs); stars: stratified permutation significance. Same-compartment pairs are excluded by design (`PairScopeFilter.CROSS_COMPARTMENT_ONLY`).
]

#let caption_supp_pan_cohort_pair_umap(pair_label, n_fdr, n_perm) = [
  *Pairwise `Harmony I2` UMAP composite for #pair_label interactions (#n_fdr verified at $"FDR"_("ACAT") < 0.10$, #n_perm additional at stratified permutation $p < 0.05$).*
  *(a)* Dual-tier cross-compartment chords and lineage cell footprints.
  *(b)* Emergent *de novo* pairwise hubs.
  *(c)* Seed-conditioned partner $Z_("ACAT")$ field.
  *(d)* Lead-pair $k$NN spotlight with $2 times 2$ biopsy response inset (Wilson $95%$ CI).
]
