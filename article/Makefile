# Article Makefile: converts SVGs to high-res PNGs and compiles Typst document

TYPST := typst
ROOT_DIR := ..
OUTPUT_DIR := $(ROOT_DIR)/output
RESULTS_DIR := $(ROOT_DIR)/results

# Use rsvg-convert (fastest, high quality) with fallback to ImageMagick
RSVG := $(shell which rsvg-convert 2>/dev/null)
MAGICK := $(shell which magick 2>/dev/null || which convert 2>/dev/null)

define convert_svg
	@if [ ! -f $(2) ] || [ $(1) -nt $(2) ]; then \
		echo "  [PNG] $(notdir $(1)) -> $(2)"; \
		if [ -n "$(RSVG)" ]; then \
			$(RSVG) -d 300 -p 300 -b white -f png $(1) -o $(2); \
		elif [ -n "$(MAGICK)" ]; then \
			$(MAGICK) -density 300 $(1) -background white -flatten $(2); \
		else \
			echo "Error: Neither rsvg-convert nor ImageMagick found for SVG conversion" && exit 1; \
		fi; \
	fi
endef

.PHONY: all article convert-figures compile-article watch clean

all: article

article: convert-figures compile-article

convert-figures:
	@echo "Ensuring figure directory structure..."
	@mkdir -p figures/dataset_overview \
	         figures/single-cell \
	         figures/deconvolution \
	         figures/mutations \
	         figures/vaf-tmb \
	         figures/immune-subsets \
	         figures/signatures \
	         figures/milo-deconv \
	         figures/coulton-tam \
	         figures/pan-cancer-atlases \
	         figures/stability-selection
	@echo "Converting dataset overview and provenance DAG SVGs to high-resolution PNGs..."
	$(call convert_svg,figures/dataset_overview/figure1_study_cohorts_overview.svg,figures/dataset_overview/figure1_study_cohorts_overview.png)
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4e_generate_tme_datasets_pipeline_dag/figure_tme_datasets_pipeline_dag.svg,figures/dataset_overview/figure_tme_datasets_pipeline_dag.png)
	@echo "Converting single-cell SVGs to high-resolution PNGs..."
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4a2_generate_combined_figures/milopy_cross_cohort_volcano_grid.svg,figures/single-cell/milopy_cross_cohort_volcano_grid.png)
	$(call convert_svg,$(OUTPUT_DIR)/single-cell-exploration/umap_cell_type_annotations.svg,figures/single-cell/umap_cell_type_annotations.png)
	$(call convert_svg,$(OUTPUT_DIR)/single-cell-exploration/qc_metadata_distributions.svg,figures/single-cell/qc_metadata_distributions.png)
	$(call convert_svg,$(OUTPUT_DIR)/single-cell-exploration/umap_granularities.svg,figures/single-cell/umap_granularities.png)
	$(call convert_svg,$(OUTPUT_DIR)/single-cell-exploration/umap_tumor_inference.svg,figures/single-cell/umap_tumor_inference.png)
	$(call convert_svg,$(OUTPUT_DIR)/single-cell-exploration/subclustering_metrics_distributions.svg,figures/single-cell/subclustering_metrics_distributions.png)
	@echo "Converting deconvolution SVGs to high-resolution PNGs..."
	$(call convert_svg,$(OUTPUT_DIR)/single-cell-exploration/deconv_prediction_advanced_grid_kmeans_subcluster_res_0.5.svg,figures/deconvolution/deconv_prediction_advanced_grid_kmeans_subcluster_res_0.5.png)
	$(call convert_svg,$(OUTPUT_DIR)/single-cell-exploration/deconv_distributions_kmeans_subcluster_res_0.5.svg,figures/deconvolution/deconv_distributions_kmeans_subcluster_res_0.5.png)
	$(call convert_svg,$(OUTPUT_DIR)/single-cell-exploration/deconv_split_response_kmeans_subcluster_res_0.5.svg,figures/deconvolution/deconv_split_response_kmeans_subcluster_res_0.5.png)
	$(call convert_svg,$(OUTPUT_DIR)/single-cell-exploration/deconv_split_tmb_kmeans_subcluster_res_0.5.svg,figures/deconvolution/deconv_split_tmb_kmeans_subcluster_res_0.5.png)
	$(call convert_svg,$(OUTPUT_DIR)/single-cell-exploration/kmeans_cluster_sizes.svg,figures/deconvolution/kmeans_cluster_sizes.png)
	$(call convert_svg,$(OUTPUT_DIR)/single-cell-exploration/kmeans_centroids_heatmap.svg,figures/deconvolution/kmeans_centroids_heatmap.png)
	$(call convert_svg,$(OUTPUT_DIR)/single-cell-exploration/iatlas_classification_umap.svg,figures/deconvolution/iatlas_classification_umap.png)
	$(call convert_svg,$(OUTPUT_DIR)/single-cell-exploration/deconv_prediction_metrics_grid.svg,figures/deconvolution/deconv_prediction_metrics_grid.png)
	@echo "Converting somatic mutation SVGs to high-resolution PNGs..."
	$(call convert_svg,$(OUTPUT_DIR)/iAtlas-response-mutations/Rosenberg-iAtlas/overfitting_underfitting_comparison.svg,figures/mutations/overfitting_underfitting_rosenberg.png)
	$(call convert_svg,$(OUTPUT_DIR)/iAtlas-response-mutations/Combined-Melanoma/overfitting_underfitting_comparison.svg,figures/mutations/overfitting_underfitting_melanoma.png)
	@echo "Converting subclonal VAF & TMB SVGs to high-resolution PNGs..."
	$(call convert_svg,$(OUTPUT_DIR)/iAtlas-vaf-tmb-stats/Rosenberg-iAtlas_vaf_tmb_association.svg,figures/vaf-tmb/Rosenberg-iAtlas_vaf_tmb_association.png)
	$(call convert_svg,$(OUTPUT_DIR)/iAtlas-vaf-tmb-reliability/vaf_density_valley_detection.svg,figures/vaf-tmb/vaf_density_valley_detection.png)
	$(call convert_svg,$(OUTPUT_DIR)/iAtlas-vaf-tmb-reliability/estimated_vs_true_thresholds.svg,figures/vaf-tmb/estimated_vs_true_thresholds.png)
	$(call convert_svg,$(OUTPUT_DIR)/iAtlas-vaf-tmb-reliability/reliability_score_vs_performance.svg,figures/vaf-tmb/reliability_score_vs_performance.png)
	$(call convert_svg,$(OUTPUT_DIR)/iAtlas-vaf-tmb-reliability/dna_vs_rna_tmb_performance.svg,figures/vaf-tmb/dna_vs_rna_tmb_performance.png)
	$(call convert_svg,$(OUTPUT_DIR)/iAtlas-vaf-tmb-optimization/odds_ratios_forest_plot_MCC.svg,figures/vaf-tmb/odds_ratios_forest_plot_MCC.png)
	$(call convert_svg,$(OUTPUT_DIR)/iAtlas-vaf-tmb-optimization/odds_ratios_forest_plot_F1.svg,figures/vaf-tmb/odds_ratios_forest_plot_F1.png)
	@echo "Copying and converting immune subset figures..."
	@cp $(OUTPUT_DIR)/results/tcell-CD8/cd8_dge_combined_figure.png figures/immune-subsets/
	@cp $(OUTPUT_DIR)/results/tcell-CD4/cd4_dge_combined_figure.png figures/immune-subsets/
	@cp $(OUTPUT_DIR)/results/myeloid/milopy_umap_response_Combined.png figures/immune-subsets/
	@cp $(OUTPUT_DIR)/results/myeloid/milopy_umap_gradient_Combined.png figures/immune-subsets/
	$(call convert_svg,$(OUTPUT_DIR)/results/myeloid/myeloid_responder_dge_volcano.svg,figures/immune-subsets/myeloid_responder_dge_volcano.png)
	$(call convert_svg,$(OUTPUT_DIR)/results/myeloid/myeloid_responder_dge_dotplot.svg,figures/immune-subsets/myeloid_responder_dge_dotplot.png)
	@echo "Converting ImmunoCompass signature SVGs to high-resolution PNGs..."
	$(call convert_svg,$(OUTPUT_DIR)/immuno-compass-signatures/signature_distributions.svg,figures/signatures/signature_distributions.png)
	$(call convert_svg,$(OUTPUT_DIR)/immuno-compass-signatures/signature_overlap_heatmap.svg,figures/signatures/signature_overlap_heatmap.png)
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4c6_generate_combined_summary_forest/fig_table_s2_combined_forest_summary.svg,figures/signatures/fig_table_s2_combined_forest_summary.png)
	@cp $(OUTPUT_DIR)/manuscript_pipeline/3e5_integrate_deconvolution_our/fig_pr_auc_heatmap_transcriptomic_pre_standard.svg figures/signatures/
	@cp $(OUTPUT_DIR)/manuscript_pipeline/3e5_integrate_deconvolution_our/fig_roc_auc_heatmap_transcriptomic_pre_standard.svg figures/signatures/
	@cp $(OUTPUT_DIR)/manuscript_pipeline/4g_plot_milo_cross_cohort_summaries/meta42_observed_recist_radar_panel_b.svg figures/milo-deconv/
	@cp $(OUTPUT_DIR)/manuscript_pipeline/4g_plot_milo_cross_cohort_summaries/pan_cohort_harmony_i2_umap_panel_b.svg figures/milo-pairwise/
	@cp $(OUTPUT_DIR)/manuscript_pipeline/4g_plot_milo_cross_cohort_summaries/milo_cross_cohort_volcano_summary.svg figures/single-cell/
	@cp $(OUTPUT_DIR)/manuscript_pipeline/4g_plot_milo_cross_cohort_summaries/milo_cross_cohort_umap_summary.svg figures/single-cell/
	@cp $(OUTPUT_DIR)/manuscript_pipeline/4i_plot_combined_main_figure/combined_main_figure_principled_gmrf.svg figures/milo-deconv/combined_main_figure.svg
	@cp $(OUTPUT_DIR)/manuscript_pipeline/4i_plot_combined_main_figure/combined_main_figure_principled_gmrf.svg figures/milo-deconv/combined_main_figure_principled_gmrf.svg
	@cp $(OUTPUT_DIR)/manuscript_pipeline/4i_plot_combined_main_figure/iatlas_milopy_vs_deconv_scatter_grid_supp_principled_gmrf.svg figures/milo-deconv/iatlas_milopy_vs_deconv_scatter_grid_supp.svg
	@cp $(OUTPUT_DIR)/manuscript_pipeline/4i_plot_combined_main_figure/iatlas_milopy_vs_deconv_scatter_grid_supp_principled_gmrf.svg figures/milo-deconv/iatlas_milopy_vs_deconv_scatter_grid_supp_principled_gmrf.svg
	@cp $(OUTPUT_DIR)/manuscript_pipeline/4j_plot_combined_markers_radar_coulton/combined_markers_radar_coulton.svg figures/milo-deconv/
	@cp $(OUTPUT_DIR)/manuscript_pipeline/4l_plot_deconv_approaches_supplementary/supp_fig_deconv_*.svg figures/milo-deconv/
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4l_plot_deconv_approaches_supplementary/supp_fig_deconv_factorial_solver_gene_matrix.svg,figures/milo-deconv/supp_fig_deconv_factorial_solver_gene_matrix.png)
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4l_plot_deconv_approaches_supplementary/supp_fig_deconv_four_axis_validation_suite.svg,figures/milo-deconv/supp_fig_deconv_four_axis_validation_suite.png)
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4l_plot_deconv_approaches_supplementary/supp_fig_deconv_approaches_head_to_head_comparison.svg,figures/milo-deconv/supp_fig_deconv_approaches_head_to_head_comparison.png)
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4l_plot_deconv_approaches_supplementary/supp_fig_deconv_gene_selection_mechanisms.svg,figures/milo-deconv/supp_fig_deconv_gene_selection_mechanisms.png)
	@echo "Converting Milo DA vs. Deconvolution SVGs to high-resolution PNGs..."
	$(call convert_svg,$(RESULTS_DIR)/sade_feldman_deconv_validation/sade_feldman_self_concordance.svg,figures/milo-deconv/sade_feldman_self_concordance.png)
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4d1_plot_sf_multires_validation_figures/fig1_predictive_capacity_benchmark.svg,figures/milo-deconv/fig1_predictive_capacity_benchmark.png)
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4d2_plot_cohort_auc_reference_selection/fig1_selected_cohort_auc_and_marginal_gain.svg,figures/milo-deconv/fig1_selected_cohort_auc_and_marginal_gain.png)
	$(call convert_svg,$(RESULTS_DIR)/synthetic_benchmark/compound_perturbation_sensitivity_curves.svg,figures/milo-deconv/compound_perturbation_sensitivity_curves.png)
	$(call convert_svg,$(RESULTS_DIR)/synthetic_benchmark/compound_perturbation_scatter_grid.svg,figures/milo-deconv/compound_perturbation_scatter_grid.png)
	$(call convert_svg,$(RESULTS_DIR)/synthetic_benchmark/factorial_2d_interaction_heatmap.svg,figures/milo-deconv/factorial_2d_interaction_heatmap.png)
	$(call convert_svg,$(RESULTS_DIR)/synthetic_benchmark/perturbation_sensitivity_curves.svg,figures/milo-deconv/perturbation_sensitivity_curves.png)
	$(call convert_svg,$(RESULTS_DIR)/synthetic_benchmark/synthetic_perturbation_scatter_grid.svg,figures/milo-deconv/synthetic_perturbation_scatter_grid.png)
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures/iatlas_response_prediction_auc_heatmap_and_bars.svg,figures/milo-deconv/iatlas_response_prediction_auc_heatmap_and_bars.png)
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures/iatlas_response_prediction_prauc_heatmap_and_bars.svg,figures/milo-deconv/iatlas_response_prediction_prauc_heatmap_and_bars.png)
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures/iatlas_milopy_vs_deconv_concordance_grid.svg,figures/milo-deconv/iatlas_milopy_vs_deconv_concordance_grid.png)
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures/iatlas_cross_cohort_concordance_transfer_matrix.svg,figures/milo-deconv/iatlas_cross_cohort_concordance_transfer_matrix.png)
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures/iatlas_consensus_metastate_forest_and_tme_profile.svg,figures/milo-deconv/iatlas_consensus_metastate_forest_and_tme_profile.png)
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures/iatlas_consensus_metastate_all42_forest_and_tme_profile.svg,figures/milo-deconv/iatlas_consensus_metastate_all42_forest_and_tme_profile.png)
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures/iatlas_consensus_metastate_top_deg_markers.svg,figures/milo-deconv/iatlas_consensus_metastate_top_deg_markers.png)
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures/iatlas_consensus_metastate_all42_top_deg_markers.svg,figures/milo-deconv/iatlas_consensus_metastate_all42_top_deg_markers.png)
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures/iatlas_classifier_and_reference_tier_pareto.svg,figures/milo-deconv/iatlas_classifier_and_reference_tier_pareto.png)
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures/iatlas_epistasis_and_stromal_gating_main.svg,figures/milo-deconv/iatlas_epistasis_and_stromal_gating_main.png)
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures/iatlas_nonlinear_conditional_interactions_all42_pure_clr.svg,figures/milo-deconv/iatlas_nonlinear_conditional_interactions_all42_pure_clr.png)
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures/iatlas_epistasis_discovery_power_benchmark.svg,figures/milo-deconv/iatlas_epistasis_discovery_power_benchmark.png)
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures/iatlas_zero_leakage_gating_comparison.svg,figures/milo-deconv/iatlas_zero_leakage_gating_comparison.png)
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures/iatlas_oof_supplementary_diagnostics.svg,figures/milo-deconv/iatlas_oof_supplementary_diagnostics.png)
	@echo "Converting Coulton TAM Atlas comparison SVGs to high-resolution PNGs..."
	@cp $(OUTPUT_DIR)/results/coulton_tam_comparison/figures/fig1b_meta42_myeloid_coulton_simplified_bicluster.svg figures/coulton-tam/
	$(call convert_svg,$(OUTPUT_DIR)/results/coulton_tam_comparison/figures/fig1_meta42_myeloid_coulton_similarity_and_ici.svg,figures/coulton-tam/fig1_meta42_myeloid_coulton_similarity_and_ici.png)
	$(call convert_svg,$(OUTPUT_DIR)/results/coulton_tam_comparison/figures/fig1b_meta42_myeloid_coulton_simplified_bicluster.svg,figures/coulton-tam/fig1b_meta42_myeloid_coulton_simplified_bicluster.png)
	$(call convert_svg,$(OUTPUT_DIR)/results/coulton_tam_comparison/figures/fig2_macrophage_nearest_substate_ici_forest.svg,figures/coulton-tam/fig2_macrophage_nearest_substate_ici_forest.png)
	$(call convert_svg,$(OUTPUT_DIR)/results/coulton_tam_comparison/figures/fig3_all42_pan_lineage_coulton_specificity.svg,figures/coulton-tam/fig3_all42_pan_lineage_coulton_specificity.png)
	@echo "Converting Pan-Cancer T-Cell, CAF & Triad Atlas comparison SVGs to high-resolution PNGs..."
	$(call convert_svg,$(OUTPUT_DIR)/figures/pan_cancer_atlas_comparison/meta42_Zheng_TCell_Atlas_similarity_and_ici.svg,figures/pan-cancer-atlases/meta42_Zheng_TCell_Atlas_similarity_and_ici.png)
	$(call convert_svg,$(OUTPUT_DIR)/figures/pan_cancer_atlas_comparison/meta42_Chu_TCell_Atlas_similarity_and_ici.svg,figures/pan-cancer-atlases/meta42_Chu_TCell_Atlas_similarity_and_ici.png)
	$(call convert_svg,$(OUTPUT_DIR)/figures/pan_cancer_atlas_comparison/meta42_Shi_PanCancer_Atlas_TNK_similarity_and_ici.svg,figures/pan-cancer-atlases/meta42_Shi_PanCancer_Atlas_TNK_similarity_and_ici.png)
	$(call convert_svg,$(OUTPUT_DIR)/figures/pan_cancer_atlas_comparison/meta42_Gao_CAF_Atlas_similarity_and_ici.svg,figures/pan-cancer-atlases/meta42_Gao_CAF_Atlas_similarity_and_ici.png)
	$(call convert_svg,$(OUTPUT_DIR)/figures/pan_cancer_atlas_comparison/meta42_Cords_CAF_Atlas_similarity_and_ici.svg,figures/pan-cancer-atlases/meta42_Cords_CAF_Atlas_similarity_and_ici.png)
	$(call convert_svg,$(OUTPUT_DIR)/figures/pan_cancer_atlas_comparison/meta42_Buechler_CAF_Atlas_similarity_and_ici.svg,figures/pan-cancer-atlases/meta42_Buechler_CAF_Atlas_similarity_and_ici.png)
	$(call convert_svg,$(OUTPUT_DIR)/figures/pan_cancer_atlas_comparison/meta42_Shi_PanCancer_Atlas_CAF_similarity_and_ici.svg,figures/pan-cancer-atlases/meta42_Shi_PanCancer_Atlas_CAF_similarity_and_ici.png)
	$(call convert_svg,$(OUTPUT_DIR)/figures/pan_cancer_atlas_comparison/meta42_triad_all42_external_specificity.svg,figures/pan-cancer-atlases/meta42_triad_all42_external_specificity.png)
	@echo "All figures successfully converted to PNG in article/figures/."

compile-article:
	@echo "Compiling main.typ into main.pdf with PNG figures..."
	$(TYPST) compile main.typ main.pdf
	@echo "Article compiled successfully: article/main.pdf"

watch:
	$(TYPST) watch main.typ main.pdf

clean:
	@rm -f main.pdf
	@rm -rf figures/*/*.png figures/*/*.svg
