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
	         figures/milo-deconv
	@echo "Converting dataset provenance DAG SVG to high-resolution PNG..."
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4d_generate_tme_datasets_pipeline_dag/figure_tme_datasets_pipeline_dag.svg,figures/dataset_overview/figure_tme_datasets_pipeline_dag.png)
	@echo "Converting single-cell SVGs to high-resolution PNGs..."
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
	@echo "Converting Milo DA vs. Deconvolution SVGs to high-resolution PNGs..."
	$(call convert_svg,$(RESULTS_DIR)/sade_feldman_deconv_validation/sade_feldman_self_concordance.svg,figures/milo-deconv/sade_feldman_self_concordance.png)
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
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures/iatlas_nonlinear_conditional_interactions.svg,figures/milo-deconv/iatlas_nonlinear_conditional_interactions.png)
	$(call convert_svg,$(OUTPUT_DIR)/manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures/iatlas_nonlinear_conditional_interactions_all42.svg,figures/milo-deconv/iatlas_nonlinear_conditional_interactions_all42.png)
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
