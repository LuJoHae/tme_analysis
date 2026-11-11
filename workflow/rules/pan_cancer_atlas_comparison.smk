# workflow/rules/pan_cancer_atlas_comparison.smk
# Snakemake rules materializing Pan-Cancer T-Cell & CAF reference atlases via tme_datasets
# and benchmarking Consensus-Meta42 across the T-Cell, CAF, and TAM Triad + 11 Bulk ICB Cohorts

import os

if "DATA_DIR" not in globals():
    DATA_DIR = os.environ.get(
        "TME_DATA_DIR", config.get("data_dir", "/storage/halu/data-test")
    )
if "RESULTS_DIR" not in globals():
    RESULTS_DIR = f"{DATA_DIR}/results"
if "MANUSCRIPT_OUT_DIR" not in globals():
    MANUSCRIPT_OUT_DIR = os.environ.get(
        "TME_PIPELINE_OUTPUT_ROOT",
        config.get("manuscript_out_dir", f"{DATA_DIR}/output/manuscript_pipeline"),
    )

PAN_ATLAS_CMP_OUT_DIR = config.get(
    "pan_atlas_cmp_out_dir", f"{RESULTS_DIR}/pan_cancer_atlas_comparison"
)
PAN_ATLAS_FIG_DIR = config.get(
    "pan_atlas_fig_dir", f"{DATA_DIR}/output/figures/pan_cancer_atlas_comparison"
)

PAN_ATLAS_FIGURE_IDS = [
    "Zheng_TCell_Atlas",
    "Chu_TCell_Atlas",
    "Shi_PanCancer_Atlas_TNK",
    "Gao_CAF_Atlas",
    "Cords_CAF_Atlas",
    "Buechler_CAF_Atlas",
    "Shi_PanCancer_Atlas_CAF",
]

PAN_ATLAS_BULK_11_COHORTS = [
    "Hugo-iAtlas",
    "Riaz-iAtlas",
    "Liu-iAtlas",
    "Gide-iAtlas",
    "Rosenberg-iAtlas",
    "Padron-iAtlas",
    "Anders-iAtlas",
    "McDermott-iAtlas",
    "Choueiri-iAtlas",
    "Jung",
    "Cho",
]

PAN_ATLAS_EXTERNAL_H5ADS = [
    "Gao_CAF_Atlas",
    "Cords_CAF_Atlas",
    "Buechler_CAF_Atlas",
    "Zheng_TCell_Atlas",
    "Shi_PanCancer_Atlas",
    "Coulton_TAM_Atlas",
]


rule all_pan_cancer_atlas_comparison:
    input:
        f"{PAN_ATLAS_CMP_OUT_DIR}/external_atlas_substate_centroids.parquet",
        f"{PAN_ATLAS_CMP_OUT_DIR}/external_atlas_top_markers.parquet",
        f"{PAN_ATLAS_CMP_OUT_DIR}/meta42_vs_external_atlases_similarity.parquet",
        f"{PAN_ATLAS_CMP_OUT_DIR}/meta42_tcell_nearest_substates.parquet",
        f"{PAN_ATLAS_CMP_OUT_DIR}/meta42_caf_nearest_substates.parquet",
        f"{PAN_ATLAS_CMP_OUT_DIR}/cross_atlas_substate_equivalence.parquet",
        f"{PAN_ATLAS_CMP_OUT_DIR}/external_atlas_substate_ici_outcomes.parquet",
        f"{PAN_ATLAS_CMP_OUT_DIR}/caf_tcell_stromal_gating_icb.parquet",
        f"{PAN_ATLAS_CMP_OUT_DIR}/triad_block_specificity_summary.parquet",
        f"{PAN_ATLAS_CMP_OUT_DIR}/external_atlas_comparison_report.md",
        expand(f"{PAN_ATLAS_FIG_DIR}/meta42_{{aid}}_similarity_and_ici.svg", aid=PAN_ATLAS_FIGURE_IDS),
        expand(f"{PAN_ATLAS_FIG_DIR}/meta42_{{aid}}_similarity_and_ici.png", aid=PAN_ATLAS_FIGURE_IDS),
        f"{PAN_ATLAS_FIG_DIR}/meta42_triad_all42_external_specificity.svg",
        f"{PAN_ATLAS_FIG_DIR}/meta42_triad_all42_external_specificity.png",


PAN_ATLAS_HUB_SUMMARY = (
    f"{MANUSCRIPT_OUT_DIR}/3c_iatlas_instaprism_milopy_benchmark/consensus_metastate_all42_hub_summary.parquet"
    if os.path.exists(
        f"{MANUSCRIPT_OUT_DIR}/3c_iatlas_instaprism_milopy_benchmark/consensus_metastate_all42_hub_summary.parquet"
    )
    else "output/results/iatlas_instaprism_benchmark/consensus_metastate_all42_hub_summary.parquet"
)


rule compare_meta42_external_atlases_tables:
    input:
        script="scripts/pan_cancer_atlas_comparison/compare_meta42_external_atlases.py",
        core_mod="scripts/pan_cancer_atlas_comparison/core.py",
        models_mod="scripts/pan_cancer_atlas_comparison/models.py",
        atlas_h5ads=expand(f"{DATA_DIR}/preprocessed/{{aid}}.h5ad", aid=PAN_ATLAS_EXTERNAL_H5ADS),
        bulk_h5ads=expand(f"{DATA_DIR}/preprocessed/{{cid}}.h5ad", cid=PAN_ATLAS_BULK_11_COHORTS),
        phi_meta=f"{MANUSCRIPT_OUT_DIR}/3b_pan_cohort_unified_deconv_benchmark/unified_reference_phi_meta.parquet",
        conc_3b=f"{MANUSCRIPT_OUT_DIR}/3b_pan_cohort_unified_deconv_benchmark/unified_all_states_detailed_concordance.parquet",
        hub_42=PAN_ATLAS_HUB_SUMMARY,
        eff_42=f"{MANUSCRIPT_OUT_DIR}/3c_iatlas_instaprism_milopy_benchmark/milopy_iatlas_state_effects.parquet",
        degs_42=f"{MANUSCRIPT_OUT_DIR}/3c_iatlas_instaprism_milopy_benchmark/consensus_metastate_all42_top_degs.parquet",
    output:
        centroids=f"{PAN_ATLAS_CMP_OUT_DIR}/external_atlas_substate_centroids.parquet",
        markers=f"{PAN_ATLAS_CMP_OUT_DIR}/external_atlas_top_markers.parquet",
        similarity=f"{PAN_ATLAS_CMP_OUT_DIR}/meta42_vs_external_atlases_similarity.parquet",
        nearest_tcell=f"{PAN_ATLAS_CMP_OUT_DIR}/meta42_tcell_nearest_substates.parquet",
        nearest_caf=f"{PAN_ATLAS_CMP_OUT_DIR}/meta42_caf_nearest_substates.parquet",
        equivalence=f"{PAN_ATLAS_CMP_OUT_DIR}/cross_atlas_substate_equivalence.parquet",
        ici_outcomes=f"{PAN_ATLAS_CMP_OUT_DIR}/external_atlas_substate_ici_outcomes.parquet",
        gating=f"{PAN_ATLAS_CMP_OUT_DIR}/caf_tcell_stromal_gating_icb.parquet",
        triad_summary=f"{PAN_ATLAS_CMP_OUT_DIR}/triad_block_specificity_summary.parquet",
        report=f"{PAN_ATLAS_CMP_OUT_DIR}/external_atlas_comparison_report.md",
    params:
        manuscript_out_dir=MANUSCRIPT_OUT_DIR,
        out_dir=PAN_ATLAS_CMP_OUT_DIR,
        fig_dir=PAN_ATLAS_FIG_DIR,
        batch_size=25000,
    threads: 4
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --manuscript-out-dir {params.manuscript_out_dir} \
            --output-dir {params.out_dir} \
            --figures-dir {params.fig_dir} \
            --batch-size {params.batch_size} \
            --skip-plots
        """


rule render_pan_cancer_atlas_comparison_figures:
    input:
        script="scripts/pan_cancer_atlas_comparison/render_figures.py",
        plotting_mod="scripts/pan_cancer_atlas_comparison/plotting.py",
        similarity=f"{PAN_ATLAS_CMP_OUT_DIR}/meta42_vs_external_atlases_similarity.parquet",
        nearest_tcell=f"{PAN_ATLAS_CMP_OUT_DIR}/meta42_tcell_nearest_substates.parquet",
        nearest_caf=f"{PAN_ATLAS_CMP_OUT_DIR}/meta42_caf_nearest_substates.parquet",
        ici_outcomes=f"{PAN_ATLAS_CMP_OUT_DIR}/external_atlas_substate_ici_outcomes.parquet",
        triad_summary=f"{PAN_ATLAS_CMP_OUT_DIR}/triad_block_specificity_summary.parquet",
    output:
        fig_atlas_svgs=expand(f"{PAN_ATLAS_FIG_DIR}/meta42_{{aid}}_similarity_and_ici.svg", aid=PAN_ATLAS_FIGURE_IDS),
        fig_atlas_pngs=expand(f"{PAN_ATLAS_FIG_DIR}/meta42_{{aid}}_similarity_and_ici.png", aid=PAN_ATLAS_FIGURE_IDS),
        fig_triad_svg=f"{PAN_ATLAS_FIG_DIR}/meta42_triad_all42_external_specificity.svg",
        fig_triad_png=f"{PAN_ATLAS_FIG_DIR}/meta42_triad_all42_external_specificity.png",
    params:
        results_dir=PAN_ATLAS_CMP_OUT_DIR,
        fig_dir=PAN_ATLAS_FIG_DIR,
    threads: 2
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --results-dir {params.results_dir} \
            --figures-dir {params.fig_dir}
        """
