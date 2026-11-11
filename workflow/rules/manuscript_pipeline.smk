# workflow/rules/manuscript_pipeline.smk
# Modular Snakemake rules orchestrating tme_datasets.pipeline + the 37-script Manuscript Pipeline (2a1-4e)
# writing exclusively to per-script isolated directories under MANUSCRIPT_OUT_DIR.

MP_OUT = os.environ.get(
    "TME_PIPELINE_OUTPUT_ROOT", config.get("manuscript_out_dir", "output/manuscript_pipeline")
)
MP_DATA_DIR = config.get("data_dir", DATA_DIR if "DATA_DIR" in locals() else "/ritom/scratch/cscs/lhaeuser/data-test")

MP_TME_INGEST_DIR = "output/tme_datasets/ingest"
MP_TME_BATCH_DIR = "output/tme_datasets/batch_preprocess"
MP_TME_QC_DIR = "output/tme_datasets/qc"

MP_2A1_DIR = f"{MP_OUT}/2a1_sample_and_analyze_cohorts"
MP_2B_DIR = f"{MP_OUT}/2b_run_stability_selection_iatlas"

MP_3A_DIR = f"{MP_OUT}/3a_cross_cohort_pseudobulk_deconv_benchmark"
MP_3B_DIR = f"{MP_OUT}/3b_pan_cohort_unified_deconv_benchmark"
MP_3C_DIR = f"{MP_OUT}/3c_iatlas_instaprism_milopy_benchmark"
MP_3D_DIR = f"{MP_OUT}/3d_curate_expert_gene_list"
MP_3E1_DIR = f"{MP_OUT}/3e1_benchmark_iatlas_icb_predictors"
MP_3E2_DIR = f"{MP_OUT}/3e2_run_hpo_pipeline"
MP_3E3_DIR = f"{MP_OUT}/3e3_benchmark_survival_and_dca"
MP_3E4_DIR = f"{MP_OUT}/3e4_update_texh_benchmark"
MP_3E5_DIR = f"{MP_OUT}/3e5_integrate_deconvolution_our"

MP_4A2_DIR = f"{MP_OUT}/4a2_generate_combined_figures"
MP_4A3_DIR = f"{MP_OUT}/4a3_plot_cross_cohort_umaps_grid"
MP_4A4_DIR = f"{MP_OUT}/4a4_plot_cross_cohort_da_bars"
MP_4A5_DIR = f"{MP_OUT}/4a5_plot_cross_cohort_deconv_scatter_grid"
MP_4A6_DIR = f"{MP_OUT}/4a6_plot_cross_cohort_deconv_correlations_bar"
MP_4A7_DIR = f"{MP_OUT}/4a7_assemble_immunotherapy_2x2_figure"
MP_4B_DIR = f"{MP_OUT}/4b_plot_iatlas_instaprism_milopy_figures"
MP_4C1_DIR = f"{MP_OUT}/4c1_plot_stability_selection_iatlas"
MP_4C2_DIR = f"{MP_OUT}/4c2_combine_all_stability_paths"
MP_4C3_DIR = f"{MP_OUT}/4c3_plot_cohort_gene_overlap"
MP_4C4_DIR = f"{MP_OUT}/4c4_plot_fitter_gene_overlap"
MP_4C5_DIR = f"{MP_OUT}/4c5_generate_hpo_summary_figures"
MP_4C6_DIR = f"{MP_OUT}/4c6_generate_combined_summary_forest"
MP_4D1_DIR = f"{MP_OUT}/4d1_plot_sf_multires_validation_figures"
MP_4D2_DIR = f"{MP_OUT}/4d2_plot_cohort_auc_reference_selection"
MP_4E_DIR = f"{MP_OUT}/4e_generate_tme_datasets_pipeline_dag"

MP_SC_9_COHORTS = [
    "GSE120575",
    "CELLxGENE_7b20c613",
    "CELLxGENE_05a8c945",
    "CELLxGENE_6f9de485",
    "GSE207422",
    "GSE243013",
    "GSE233203",
    "GSE200996",
    "GSE316195",
]

MP_STABILITY_12_SLUGS = [
    "pancancer",
    "melanoma",
    "rcc",
    "rosenberg_iatlas",
    "mcdermott_iatlas",
    "riaz_iatlas",
    "gide_iatlas",
    "liu_iatlas",
    "hugo_iatlas",
    "choueiri_iatlas",
    "padron_iatlas",
    "anders_iatlas",
]


# ==============================================================================
# Aggregate Target: Complete 37-Script Manuscript Pipeline (2a1-4e)
# ==============================================================================
rule all_manuscript_pipeline:
    input:
        # Tier IV [4A]: Single-cell 9-cohort volcano grid, UMAPs, pseudobulk scatters & 2x2 master
        f"{MP_4A2_DIR}/milopy_cross_cohort_volcano_grid.svg",
        f"{MP_4A2_DIR}/milopy_cross_cohort_volcano_grid.png",
        f"{MP_4A3_DIR}/nine_cohorts_umaps_logfc_grid.svg",
        f"{MP_4A3_DIR}/nine_cohorts_umaps_logfc_grid.png",
        f"{MP_4A4_DIR}/nine_cohorts_da_bars.svg",
        f"{MP_4A4_DIR}/nine_cohorts_da_bars.png",
        f"{MP_4A5_DIR}/nine_cohorts_deconv_scatter_grid.svg",
        f"{MP_4A5_DIR}/nine_cohorts_deconv_scatter_grid.png",
        f"{MP_4A6_DIR}/nine_cohorts_deconv_correlations_bar.svg",
        f"{MP_4A6_DIR}/nine_cohorts_deconv_composite_with_bars.svg",
        f"{MP_4A7_DIR}/immunotherapy_response_2x2_master.svg",
        f"{MP_4A7_DIR}/immunotherapy_response_2x2_master.png",
        # Tier IV [4B]: 14 iAtlas InstaPrism x Milopy concordance, Main Epistasis/Gating & Supp S19-S22 figures + Typst captions
        f"{MP_4B_DIR}/iatlas_milopy_vs_deconv_concordance_grid.svg",
        f"{MP_4B_DIR}/iatlas_milopy_vs_deconv_concordance_grid.png",
        f"{MP_4B_DIR}/iatlas_response_prediction_auc_heatmap_and_bars.svg",
        f"{MP_4B_DIR}/iatlas_response_prediction_prauc_heatmap_and_bars.svg",
        f"{MP_4B_DIR}/iatlas_classifier_and_reference_tier_pareto.svg",
        f"{MP_4B_DIR}/iatlas_consensus_metastate_all42_forest_and_tme_profile.svg",
        f"{MP_4B_DIR}/iatlas_consensus_metastate_all42_top_deg_markers.svg",
        f"{MP_4B_DIR}/iatlas_epistasis_and_stromal_gating_main.svg",
        f"{MP_4B_DIR}/iatlas_epistasis_and_stromal_gating_main.png",
        f"{MP_4B_DIR}/iatlas_nonlinear_conditional_interactions_all42_pure_clr.svg",
        f"{MP_4B_DIR}/iatlas_epistasis_discovery_power_benchmark.svg",
        f"{MP_4B_DIR}/iatlas_zero_leakage_gating_comparison.svg",
        f"{MP_4B_DIR}/iatlas_oof_supplementary_diagnostics.svg",
        f"{MP_4B_DIR}/iatlas_epistasis_and_stromal_gating_main_caption.typ",
        f"{MP_4B_DIR}/iatlas_epistasis_supplementary_captions.typ",
        # Tier III/IV [3E5 & 4C]: Stability selection, Optuna HPO, 28-RNA PR-AUC heatmap & 32-biomarker forest
        f"{MP_4C2_DIR}/all_cohorts_stability_paths_composite.svg",
        f"{MP_4C2_DIR}/all_cohorts_stability_paths_composite.png",
        f"{MP_4C3_DIR}/cohort_gene_overlap_nature.svg",
        f"{MP_4C3_DIR}/cohort_gene_overlap_nature.png",
        f"{MP_4C4_DIR}/fitter_gene_overlap_nature.svg",
        f"{MP_4C4_DIR}/fitter_gene_overlap_nature.png",
        f"{MP_4C5_DIR}/fig_hpo_method_and_sampling.svg",
        f"{MP_4C5_DIR}/fig_hpo_benchmark_results.svg",
        f"{MP_3E5_DIR}/fig_pr_auc_heatmap_transcriptomic_pre_standard.svg",
        f"{MP_3E5_DIR}/fig_pr_auc_heatmap_transcriptomic_pre_standard.png",
        f"{MP_4C6_DIR}/fig_table_s2_combined_forest_summary.svg",
        f"{MP_4C6_DIR}/fig_table_s2_combined_forest_summary.png",
        # Tier IV [4D]: Multi-resolution saturation grid (a)-(r) & selected cohort AUC gain
        f"{MP_4D1_DIR}/fig1_predictive_capacity_benchmark.svg",
        f"{MP_4D1_DIR}/fig1_predictive_capacity_benchmark.png",
        f"{MP_4D2_DIR}/fig1_selected_cohort_auc_and_marginal_gain.svg",
        f"{MP_4D2_DIR}/fig1_selected_cohort_auc_and_marginal_gain.png",
        # Tier IV [4E]: 37-script Provenance DAG
        f"{MP_4E_DIR}/figure_tme_datasets_pipeline_dag.svg",
        f"{MP_4E_DIR}/figure_tme_datasets_pipeline_dag.png",
        # Curated gene dossier
        f"{MP_3D_DIR}/curated_expert_gene_dossier.parquet",


# ==============================================================================
# Tier I (packages/tme_datasets): Ingestion, Batch Preprocessing & Single-Cell QC
# ==============================================================================
rule mp_tme_datasets_ingest:
    input:
        script="packages/tme_datasets/src/tme_datasets/pipeline/download_preprocess.py",
    output:
        progress=f"{MP_TME_INGEST_DIR}/preprocessing_progress.json",
    threads: 8
    shell:
        """
        set -euo pipefail
        if [ "${{TME_DISABLE_DOWNLOADS:-0}}" = "1" ] || [ -d /ritom/scratch ]; then
            echo "ERROR: Internet dataset downloads are strictly forbidden on CSCS Eiger. Download on olm and transfer via Globus." >&2
            exit 1
        fi
        python -m tme_datasets.pipeline.download_preprocess --workers {threads}
        """


rule mp_tme_datasets_batch_preprocess:
    input:
        script="packages/tme_datasets/src/tme_datasets/pipeline/batch_preprocess.py",
    output:
        sentinel=f"{MP_TME_BATCH_DIR}/batch_preprocess_completed.txt",
    params:
        registry="data/registry/discovered_solid_tumor_sc_datasets.parquet",
        raw_dir=f"{MP_DATA_DIR}/raw",
        preprocessed_dir=f"{MP_DATA_DIR}/preprocessed",
        report_dir=MP_TME_BATCH_DIR,
    shell:
        """
        set -euo pipefail
        mkdir -p {params.report_dir}
        python -m tme_datasets.pipeline.batch_preprocess \
            --registry {params.registry} \
            --data-raw-dir {params.raw_dir} \
            --output-dir {params.preprocessed_dir} \
            --report-dir {params.report_dir}
        touch {output.sentinel}
        """


rule mp_tme_datasets_qc:
    input:
        script="packages/tme_datasets/src/tme_datasets/pipeline/qc_orchestrator.py",
    output:
        summary_csv=f"{MP_TME_QC_DIR}/single_cell_qc_run_summary.csv",
        summary_md=f"{MP_TME_QC_DIR}/single_cell_qc_run_summary.md",
    params:
        out_dir=MP_TME_QC_DIR,
    shell:
        """
        set -euo pipefail
        python -m tme_datasets.pipeline.qc_orchestrator --output-dir {params.out_dir}
        """


# ==============================================================================
# Tier III: Pseudobulk Deconvolution, Consensus-Meta42 & iAtlas Concordance (3a-3d)
# ==============================================================================
rule mp_3a_cross_cohort_pseudobulk_deconv_benchmark:
    input:
        script="scripts/manuscript_pipeline/3a_cross_cohort_pseudobulk_deconv_benchmark.py",
        da_parquets=expand(f"{MP_2A1_DIR}/{{sc_cohort}}/da_results.parquet", sc_cohort=MP_SC_9_COHORTS),
    output:
        summary=f"{MP_3A_DIR}/cross_cohort_concordance_summary.parquet",
        sweeps=f"{MP_3A_DIR}/cross_cohort_baseline_sweeps_summary.parquet",
        svg=f"{MP_3A_DIR}/cross_cohort_pseudobulk_deconv_summary.svg",
        png=f"{MP_3A_DIR}/cross_cohort_pseudobulk_deconv_summary.png",
    params:
        milo_dir=MP_2A1_DIR,
        out_dir=MP_3A_DIR,
        adata_dir=f"{MP_DATA_DIR}/preprocessed",
        pb_cache_dir=f"{MP_DATA_DIR}/pseudobulk",
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --cohort all \
            --adata-dir {params.adata_dir} \
            --pseudobulk-cache-dir {params.pb_cache_dir} \
            --milo-dir {params.milo_dir} \
            --out-dir {params.out_dir} \
            --results-dir {params.out_dir}
        """


rule mp_3b_pan_cohort_unified_deconv_benchmark:
    input:
        script="scripts/manuscript_pipeline/3b_pan_cohort_unified_deconv_benchmark.py",
        stage3a_summary=f"{MP_3A_DIR}/cross_cohort_concordance_summary.parquet",
    output:
        phi_meta=f"{MP_3B_DIR}/unified_reference_phi_meta.parquet",
        meta_conc=f"{MP_3B_DIR}/unified_meta_state_concordance.parquet",
        cohort_conc=f"{MP_3B_DIR}/unified_per_cohort_concordance_summary.parquet",
        detailed=f"{MP_3B_DIR}/unified_all_states_detailed_concordance.parquet",
        svg=f"{MP_3B_DIR}/pan_cohort_unified_reference_comparison.svg",
        png=f"{MP_3B_DIR}/pan_cohort_unified_reference_comparison.png",
    params:
        bench_dir=MP_3A_DIR,
        pb_cache_dir=f"{MP_DATA_DIR}/pseudobulk",
        out_dir=MP_3B_DIR,
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --benchmark-dir {params.bench_dir} \
            --pseudobulk-cache-dir {params.pb_cache_dir} \
            --out-dir {params.out_dir} \
            --results-dir {params.out_dir}
        """


rule mp_3c_iatlas_instaprism_milopy_benchmark:
    input:
        script="scripts/manuscript_pipeline/3c_iatlas_instaprism_milopy_benchmark.py",
        stage3a_summary=f"{MP_3A_DIR}/cross_cohort_concordance_summary.parquet",
        stage3b_phi=f"{MP_3B_DIR}/unified_reference_phi_meta.parquet",
    output:
        state_effects=f"{MP_3C_DIR}/milopy_iatlas_state_effects.parquet",
        conc_summary=f"{MP_3C_DIR}/milopy_iatlas_concordance_summary.parquet",
        resp_metrics=f"{MP_3C_DIR}/response_prediction_metrics.parquet",
        loco_metrics=f"{MP_3C_DIR}/loco_prediction_metrics.parquet",
        degs14=f"{MP_3C_DIR}/consensus_metastate_top_degs.parquet",
        degs42=f"{MP_3C_DIR}/consensus_metastate_all42_top_degs.parquet",
        inter42_pure=f"{MP_3C_DIR}/consensus_metastate_all42_interactions_pure_clr.parquet",
        zero_leak_summary=f"{MP_3C_DIR}/zero_leakage_gating_methods_summary.parquet",
        zero_leak_cells=f"{MP_3C_DIR}/zero_leakage_gating_regimes_cells.parquet",
        oof_2x2_cells=f"{MP_3C_DIR}/oof_2x2_extreme_regimes_cells.parquet",
        oof_2x2_contrasts=f"{MP_3C_DIR}/oof_2x2_extreme_regimes_contrasts.parquet",
        epi_strategies=f"{MP_3C_DIR}/epistasis_discovery_strategies_summary.parquet",
        epi_sig_pairs=f"{MP_3C_DIR}/epistasis_significant_pairs_calibrated.parquet",
        oof_9a9b_landscape=f"{MP_3C_DIR}/oof_9a_9b_patient_landscape_and_curves.parquet",
        oof_9a9b_clin=f"{MP_3C_DIR}/oof_9a_9b_clinical_and_histology_breakdown.parquet",
        oof_9a9b_drivers=f"{MP_3C_DIR}/oof_9a_9b_fold_state_drivers.parquet",
    params:
        sc_bench_dir=MP_3A_DIR,
        preprocessed_dir=f"{MP_DATA_DIR}/preprocessed",
        out_dir=MP_3C_DIR,
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --sc-benchmark-dir {params.sc_bench_dir} \
            --preprocessed-dir {params.preprocessed_dir} \
            --out-dir {params.out_dir} \
            --results-dir {params.out_dir}
        """


rule mp_3d_curate_expert_gene_list:
    input:
        script="scripts/manuscript_pipeline/3d_curate_expert_gene_list.py",
        scores=f"{MP_2B_DIR}/stability_scores.parquet",
    output:
        parquet=f"{MP_3D_DIR}/curated_expert_gene_dossier.parquet",
        md=f"{MP_3D_DIR}/expert_gene_dossier.md",
        csv=f"{MP_3D_DIR}/expert_gene_summary_table.csv",
    shell:
        """
        set -euo pipefail
        python {input.script}
        """


# ==============================================================================
# Tier III: 32-Biomarker ICB Response, Survival/DCA & Deconvolution Integration (3e1, 3e3-3e5)
# ==============================================================================
rule mp_3e1_benchmark_iatlas_icb_predictors:
    input:
        script="scripts/manuscript_pipeline/3e1_benchmark_iatlas_icb_predictors.py",
    output:
        metrics=f"{MP_3E1_DIR}/iatlas_benchmark_metrics.parquet",
        predictability=f"{MP_3E1_DIR}/cohort_predictability_metrics.parquet",
        perturbation=f"{MP_3E1_DIR}/perturbation_benchmark_metrics.parquet",
    params:
        out_dir=MP_3E1_DIR,
    shell:
        """
        set -euo pipefail
        python {input.script} --output-dir {params.out_dir} --enable-combined-pools --run-perturbations
        """


rule mp_3e3_benchmark_survival_and_dca:
    input:
        script="scripts/manuscript_pipeline/3e3_benchmark_survival_and_dca.py",
    output:
        survival=f"{MP_3E3_DIR}/survival_benchmark_metrics.parquet",
        dca=f"{MP_3E3_DIR}/dca_benchmark_metrics.parquet",
    params:
        out_dir=MP_3E3_DIR,
    shell:
        """
        set -euo pipefail
        python {input.script} --output-dir {params.out_dir} --enable-combined-pools
        """


rule mp_3e4_update_texh_benchmark:
    input:
        script="scripts/manuscript_pipeline/3e4_update_texh_benchmark.py",
        metrics=f"{MP_3E1_DIR}/iatlas_benchmark_metrics.parquet",
        survival=f"{MP_3E3_DIR}/survival_benchmark_metrics.parquet",
        dca=f"{MP_3E3_DIR}/dca_benchmark_metrics.parquet",
    output:
        metrics=f"{MP_3E4_DIR}/iatlas_benchmark_metrics.parquet",
        survival=f"{MP_3E4_DIR}/survival_benchmark_metrics.parquet",
        dca=f"{MP_3E4_DIR}/dca_benchmark_metrics.parquet",
    params:
        upstream_dir=MP_3E1_DIR,
        out_dir=MP_3E4_DIR,
    shell:
        """
        set -euo pipefail
        python {input.script} --upstream-dir {params.upstream_dir} --output-dir {params.out_dir}
        """


rule mp_3e5_integrate_deconvolution_our:
    input:
        script="scripts/manuscript_pipeline/3e5_integrate_deconvolution_our.py",
        metrics=f"{MP_3E4_DIR}/iatlas_benchmark_metrics.parquet",
        survival=f"{MP_3E4_DIR}/survival_benchmark_metrics.parquet",
        dca=f"{MP_3E4_DIR}/dca_benchmark_metrics.parquet",
        predictability=f"{MP_3E1_DIR}/cohort_predictability_metrics.parquet",
        perturbation=f"{MP_3E1_DIR}/perturbation_benchmark_metrics.parquet",
    output:
        metrics=f"{MP_3E5_DIR}/iatlas_benchmark_metrics.parquet",
        survival=f"{MP_3E5_DIR}/survival_benchmark_metrics.parquet",
        dca=f"{MP_3E5_DIR}/dca_benchmark_metrics.parquet",
        pr_rna_svg=f"{MP_3E5_DIR}/fig_pr_auc_heatmap_transcriptomic_pre_standard.svg",
        pr_rna_png=f"{MP_3E5_DIR}/fig_pr_auc_heatmap_transcriptomic_pre_standard.png",
        roc_rna_svg=f"{MP_3E5_DIR}/fig_roc_auc_heatmap_transcriptomic_pre_standard.svg",
        roc_rna_png=f"{MP_3E5_DIR}/fig_roc_auc_heatmap_transcriptomic_pre_standard.png",
        forest_svg=f"{MP_3E5_DIR}/fig_table_s2_combined_forest_summary.svg",
        forest_png=f"{MP_3E5_DIR}/fig_table_s2_combined_forest_summary.png",
    params:
        upstream_dir=MP_3E4_DIR,
        out_dir=MP_3E5_DIR,
    shell:
        """
        set -euo pipefail
        python {input.script} --upstream-dir {params.upstream_dir} --output-dir {params.out_dir}
        """


# ==============================================================================
# Tier IV: Single-Cell 9-Cohort Grids & 2x2 Master Figure (4a3-4a7)
# ==============================================================================
rule mp_4a3_plot_cross_cohort_umaps_grid:
    input:
        script="scripts/manuscript_pipeline/4a3_plot_cross_cohort_umaps_grid.py",
        cell_parquets=expand(f"{MP_2A1_DIR}/{{sc_cohort}}/cell_level_scores.parquet", sc_cohort=MP_SC_9_COHORTS),
    output:
        svg=f"{MP_4A3_DIR}/nine_cohorts_umaps_logfc_grid.svg",
        png=f"{MP_4A3_DIR}/nine_cohorts_umaps_logfc_grid.png",
    params:
        results_dir=MP_2A1_DIR,
    shell:
        """
        set -euo pipefail
        python {input.script} --results-dir {params.results_dir} --out-svg {output.svg} --out-png {output.png}
        """


rule mp_4a4_plot_cross_cohort_da_bars:
    input:
        script="scripts/manuscript_pipeline/4a4_plot_cross_cohort_da_bars.py",
        da_parquets=expand(f"{MP_2A1_DIR}/{{sc_cohort}}/da_results.parquet", sc_cohort=MP_SC_9_COHORTS),
    output:
        svg=f"{MP_4A4_DIR}/nine_cohorts_da_bars.svg",
        png=f"{MP_4A4_DIR}/nine_cohorts_da_bars.png",
    params:
        results_dir=MP_2A1_DIR,
    shell:
        """
        set -euo pipefail
        python {input.script} --results-dir {params.results_dir} --out-svg {output.svg} --out-png {output.png}
        """


rule mp_4a5_plot_cross_cohort_deconv_scatter_grid:
    input:
        script="scripts/manuscript_pipeline/4a5_plot_cross_cohort_deconv_scatter_grid.py",
        summary=f"{MP_3A_DIR}/cross_cohort_concordance_summary.parquet",
    output:
        svg=f"{MP_4A5_DIR}/nine_cohorts_deconv_scatter_grid.svg",
        png=f"{MP_4A5_DIR}/nine_cohorts_deconv_scatter_grid.png",
    params:
        bench_dir=MP_3A_DIR,
    shell:
        """
        set -euo pipefail
        python {input.script} --benchmark-dir {params.bench_dir} --out-svg {output.svg} --out-png {output.png}
        """


rule mp_4a6_plot_cross_cohort_deconv_correlations_bar:
    input:
        script="scripts/manuscript_pipeline/4a6_plot_cross_cohort_deconv_correlations_bar.py",
        summary=f"{MP_3A_DIR}/cross_cohort_concordance_summary.parquet",
        scatter_svg=f"{MP_4A5_DIR}/nine_cohorts_deconv_scatter_grid.svg",
    output:
        bar_svg=f"{MP_4A6_DIR}/nine_cohorts_deconv_correlations_bar.svg",
        bar_png=f"{MP_4A6_DIR}/nine_cohorts_deconv_correlations_bar.png",
        comp_svg=f"{MP_4A6_DIR}/nine_cohorts_deconv_composite_with_bars.svg",
        comp_png=f"{MP_4A6_DIR}/nine_cohorts_deconv_composite_with_bars.png",
    params:
        bench_dir=MP_3A_DIR,
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --benchmark-dir {params.bench_dir} \
            --scatter-svg {input.scatter_svg} \
            --out-bar-svg {output.bar_svg} \
            --out-bar-png {output.bar_png} \
            --out-composite-svg {output.comp_svg} \
            --out-composite-png {output.comp_png}
        """


rule mp_4a7_assemble_immunotherapy_2x2_figure:
    input:
        script="scripts/manuscript_pipeline/4a7_assemble_immunotherapy_2x2_figure.py",
        panel_a=f"{MP_4A3_DIR}/nine_cohorts_umaps_logfc_grid.svg",
        panel_b=f"{MP_4A4_DIR}/nine_cohorts_da_bars.svg",
        panel_c=f"{MP_4A5_DIR}/nine_cohorts_deconv_scatter_grid.svg",
        panel_d=f"{MP_4A6_DIR}/nine_cohorts_deconv_correlations_bar.svg",
    output:
        svg=f"{MP_4A7_DIR}/immunotherapy_response_2x2_master.svg",
        png=f"{MP_4A7_DIR}/immunotherapy_response_2x2_master.png",
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --panel-a {input.panel_a} \
            --panel-b {input.panel_b} \
            --panel-c {input.panel_c} \
            --panel-d {input.panel_d} \
            --out-svg {output.svg} \
            --out-png {output.png}
        """


# ==============================================================================
# Tier IV: iAtlas InstaPrism x Milopy, Stability Selection, HPO & Forest Figures (4b, 4c1-4c6)
# ==============================================================================
rule mp_4b_plot_iatlas_instaprism_milopy_figures:
    input:
        script="scripts/manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures.py",
        state_effects=f"{MP_3C_DIR}/milopy_iatlas_state_effects.parquet",
        conc_summary=f"{MP_3C_DIR}/milopy_iatlas_concordance_summary.parquet",
        resp_metrics=f"{MP_3C_DIR}/response_prediction_metrics.parquet",
        loco_metrics=f"{MP_3C_DIR}/loco_prediction_metrics.parquet",
        degs14=f"{MP_3C_DIR}/consensus_metastate_top_degs.parquet",
        degs42=f"{MP_3C_DIR}/consensus_metastate_all42_top_degs.parquet",
        inter42_pure=f"{MP_3C_DIR}/consensus_metastate_all42_interactions_pure_clr.parquet",
        zero_leak_summary=f"{MP_3C_DIR}/zero_leakage_gating_methods_summary.parquet",
        zero_leak_cells=f"{MP_3C_DIR}/zero_leakage_gating_regimes_cells.parquet",
        oof_2x2_cells=f"{MP_3C_DIR}/oof_2x2_extreme_regimes_cells.parquet",
        oof_2x2_contrasts=f"{MP_3C_DIR}/oof_2x2_extreme_regimes_contrasts.parquet",
        epi_strategies=f"{MP_3C_DIR}/epistasis_discovery_strategies_summary.parquet",
        epi_sig_pairs=f"{MP_3C_DIR}/epistasis_significant_pairs_calibrated.parquet",
        oof_9a9b_landscape=f"{MP_3C_DIR}/oof_9a_9b_patient_landscape_and_curves.parquet",
        oof_9a9b_clin=f"{MP_3C_DIR}/oof_9a_9b_clinical_and_histology_breakdown.parquet",
        oof_9a9b_drivers=f"{MP_3C_DIR}/oof_9a_9b_fold_state_drivers.parquet",
    output:
        grid_svg=f"{MP_4B_DIR}/iatlas_milopy_vs_deconv_concordance_grid.svg",
        grid_png=f"{MP_4B_DIR}/iatlas_milopy_vs_deconv_concordance_grid.png",
        matrix_svg=f"{MP_4B_DIR}/iatlas_cross_cohort_concordance_transfer_matrix.svg",
        auc_svg=f"{MP_4B_DIR}/iatlas_response_prediction_auc_heatmap_and_bars.svg",
        prauc_svg=f"{MP_4B_DIR}/iatlas_response_prediction_prauc_heatmap_and_bars.svg",
        pareto_svg=f"{MP_4B_DIR}/iatlas_classifier_and_reference_tier_pareto.svg",
        forest42_svg=f"{MP_4B_DIR}/iatlas_consensus_metastate_all42_forest_and_tme_profile.svg",
        deg42_svg=f"{MP_4B_DIR}/iatlas_consensus_metastate_all42_top_deg_markers.svg",
        main_epi_svg=f"{MP_4B_DIR}/iatlas_epistasis_and_stromal_gating_main.svg",
        main_epi_png=f"{MP_4B_DIR}/iatlas_epistasis_and_stromal_gating_main.png",
        inter42_pure_svg=f"{MP_4B_DIR}/iatlas_nonlinear_conditional_interactions_all42_pure_clr.svg",
        epi_power_svg=f"{MP_4B_DIR}/iatlas_epistasis_discovery_power_benchmark.svg",
        zero_leak_svg=f"{MP_4B_DIR}/iatlas_zero_leakage_gating_comparison.svg",
        oof_diag_svg=f"{MP_4B_DIR}/iatlas_oof_supplementary_diagnostics.svg",
        main_caption_typ=f"{MP_4B_DIR}/iatlas_epistasis_and_stromal_gating_main_caption.typ",
        supp_captions_typ=f"{MP_4B_DIR}/iatlas_epistasis_supplementary_captions.typ",
    params:
        data_dir=MP_3C_DIR,
        results_dir=MP_4B_DIR,
    shell:
        """
        set -euo pipefail
        python {input.script} --data-dir {params.data_dir} --results-dir {params.results_dir}
        """


rule mp_4c1_plot_stability_selection_iatlas:
    input:
        script="scripts/manuscript_pipeline/4c1_plot_stability_selection_iatlas.py",
        scores=f"{MP_2B_DIR}/stability_scores.parquet",
        paths=f"{MP_2B_DIR}/stability_paths.parquet",
    output:
        svgs=expand(f"{MP_4C1_DIR}/paths_comparison_{{slug}}.svg", slug=MP_STABILITY_12_SLUGS),
    params:
        in_dir=MP_2B_DIR,
        out_dir=MP_4C1_DIR,
    shell:
        """
        set -euo pipefail
        python {input.script} --input-dir {params.in_dir} --output-dir {params.out_dir}
        """


rule mp_4c2_combine_all_stability_paths:
    input:
        script="scripts/manuscript_pipeline/4c2_combine_all_stability_paths.py",
        svgs=expand(f"{MP_4C1_DIR}/paths_comparison_{{slug}}.svg", slug=MP_STABILITY_12_SLUGS),
    output:
        svg=f"{MP_4C2_DIR}/all_cohorts_stability_paths_composite.svg",
        png=f"{MP_4C2_DIR}/all_cohorts_stability_paths_composite.png",
    params:
        in_dir=MP_4C1_DIR,
    shell:
        """
        set -euo pipefail
        python {input.script} --input-dir {params.in_dir} --output-svg {output.svg} --output-png {output.png}
        """


rule mp_4c3_plot_cohort_gene_overlap:
    input:
        script="scripts/manuscript_pipeline/4c3_plot_cohort_gene_overlap.py",
        scores=f"{MP_2B_DIR}/stability_scores.parquet",
        summary=f"{MP_2B_DIR}/cohort_fitter_summary.parquet",
    output:
        svg=f"{MP_4C3_DIR}/cohort_gene_overlap_nature.svg",
        png=f"{MP_4C3_DIR}/cohort_gene_overlap_nature.png",
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --scores-path {input.scores} \
            --summary-path {input.summary} \
            --output-svg {output.svg} \
            --output-png {output.png}
        """


rule mp_4c5_generate_hpo_summary_figures:
    input:
        script="scripts/manuscript_pipeline/4c5_generate_hpo_summary_figures.py",
        hpo_evals=f"{MP_3E2_DIR}/hpo_evaluations.parquet",
        bench_metrics=f"{MP_3E1_DIR}/iatlas_benchmark_metrics.parquet",
    output:
        method_svg=f"{MP_4C5_DIR}/fig_hpo_method_and_sampling.svg",
        method_png=f"{MP_4C5_DIR}/fig_hpo_method_and_sampling.png",
        results_svg=f"{MP_4C5_DIR}/fig_hpo_benchmark_results.svg",
        results_png=f"{MP_4C5_DIR}/fig_hpo_benchmark_results.png",
    shell:
        """
        set -euo pipefail
        python {input.script}
        """


rule mp_4c6_generate_combined_summary_forest:
    input:
        script="scripts/manuscript_pipeline/4c6_generate_combined_summary_forest.py",
        metrics=f"{MP_3E5_DIR}/iatlas_benchmark_metrics.parquet",
    output:
        svg=f"{MP_4C6_DIR}/fig_table_s2_combined_forest_summary.svg",
        png=f"{MP_4C6_DIR}/fig_table_s2_combined_forest_summary.png",
    params:
        out_dir=MP_4C6_DIR,
    shell:
        """
        set -euo pipefail
        python {input.script} --parquet {input.metrics} --output-dir {params.out_dir}
        """


# ==============================================================================
# Tier IV Terminal Stage: 37-Script Provenance DAG (4e)
# ==============================================================================
rule mp_4e_generate_tme_datasets_pipeline_dag:
    input:
        script="scripts/manuscript_pipeline/4e_generate_tme_datasets_pipeline_dag.py",
    output:
        svg=f"{MP_4E_DIR}/figure_tme_datasets_pipeline_dag.svg",
        png=f"{MP_4E_DIR}/figure_tme_datasets_pipeline_dag.png",
    shell:
        """
        set -euo pipefail
        python {input.script}
        """


# ==============================================================================
# Pairwise Milopy Differential Neighborhood Co-Association Across 9 Cohorts
# ==============================================================================
MP_PAIRWISE_OUT_DIR = f"{MP_OUT}/milopy_pairwise_analysis"

rule mp_analyze_pairwise_cohorts:
    input:
        script="scripts/milopy_analysis/analyze_pairwise_cohorts.py",
        da_parquets=expand(f"{MP_2A1_DIR}/{{sc_cohort}}/da_results.parquet", sc_cohort=MP_SC_9_COHORTS),
    output:
        summary=f"{MP_PAIRWISE_OUT_DIR}/reports/pairwise_cohorts_summary.parquet",
        recurrent=f"{MP_PAIRWISE_OUT_DIR}/reports/recurrent_celltype_pairs.parquet",
        report_md=f"{MP_PAIRWISE_OUT_DIR}/reports/milopy_pairwise_cohorts_summary.md",
        comp_svg=f"{MP_PAIRWISE_OUT_DIR}/reports/pairwise_cohort_comparison.svg",
        inter_svg=f"{MP_PAIRWISE_OUT_DIR}/reports/pairwise_celltype_interactions.svg",
    params:
        in_dir=MP_2A1_DIR,
        alt_in_dir=f"{MP_DATA_DIR}/results/milopy_response_analysis",
        out_dir=MP_PAIRWISE_OUT_DIR,
    threads: 8
    shell:
        """
        set -euo pipefail
        INPUT_DIR="{params.in_dir}"
        if [ ! -f "{params.in_dir}/GSE120575/nhood_matrix.npz" ] && [ -f "{params.alt_in_dir}/GSE120575/nhood_matrix.npz" ]; then
            INPUT_DIR="{params.alt_in_dir}"
        fi
        python {input.script} --input-dir "$INPUT_DIR" --output-dir {params.out_dir}
        """
