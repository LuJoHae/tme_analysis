# workflow/rules/coulton_tam_comparison.smk
# Snakemake rules comparing Consensus-Meta42 states against Coulton_TAM_Atlas substates & ICI outcomes

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

COULTON_CMP_OUT_DIR = config.get(
    "coulton_cmp_out_dir", f"{RESULTS_DIR}/coulton_tam_comparison"
)

COULTON_BULK_11_COHORTS = [
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


rule all_coulton_tam_comparison:
    input:
        f"{COULTON_CMP_OUT_DIR}/coulton_substate_centroids.parquet",
        f"{COULTON_CMP_OUT_DIR}/coulton_top_markers.parquet",
        f"{COULTON_CMP_OUT_DIR}/meta42_vs_coulton_full_similarity.parquet",
        f"{COULTON_CMP_OUT_DIR}/meta42_macrophage_nearest_coulton_substates.parquet",
        f"{COULTON_CMP_OUT_DIR}/meta42_coulton_sensitivity_summary.parquet",
        f"{COULTON_CMP_OUT_DIR}/coulton_substate_ici_outcomes.parquet",
        f"{COULTON_CMP_OUT_DIR}/meta42_macrophage_indication_effects.parquet",
        f"{COULTON_CMP_OUT_DIR}/meta42_coulton_tam_comparison_report.md",
        f"{COULTON_CMP_OUT_DIR}/figures/fig1_meta42_myeloid_coulton_similarity_and_ici.svg",
        f"{COULTON_CMP_OUT_DIR}/figures/fig1_meta42_myeloid_coulton_similarity_and_ici.png",
        f"{COULTON_CMP_OUT_DIR}/figures/fig1b_meta42_myeloid_coulton_simplified_bicluster.svg",
        f"{COULTON_CMP_OUT_DIR}/figures/fig1b_meta42_myeloid_coulton_simplified_bicluster.png",
        f"{COULTON_CMP_OUT_DIR}/figures/fig2_macrophage_nearest_substate_ici_forest.svg",
        f"{COULTON_CMP_OUT_DIR}/figures/fig2_macrophage_nearest_substate_ici_forest.png",
        f"{COULTON_CMP_OUT_DIR}/figures/fig3_all42_pan_lineage_coulton_specificity.svg",
        f"{COULTON_CMP_OUT_DIR}/figures/fig3_all42_pan_lineage_coulton_specificity.png",
        f"{COULTON_CMP_OUT_DIR}/coulton_tam_simplified_bicluster_main_caption.typ",


rule compare_meta42_coulton_tam_tables:
    input:
        script="scripts/coulton_tam_comparison/compare_meta42_coulton_tam.py",
        core_mod="scripts/coulton_tam_comparison/core.py",
        models_mod="scripts/coulton_tam_comparison/models.py",
        coulton_h5ad=f"{DATA_DIR}/preprocessed/Coulton_TAM_Atlas.h5ad",
        bulk_h5ads=expand(f"{DATA_DIR}/preprocessed/{{cid}}.h5ad", cid=COULTON_BULK_11_COHORTS),
        phi_meta=f"{MANUSCRIPT_OUT_DIR}/3b_pan_cohort_unified_deconv_benchmark/unified_reference_phi_meta.parquet",
        conc_3b=f"{MANUSCRIPT_OUT_DIR}/3b_pan_cohort_unified_deconv_benchmark/unified_all_states_detailed_concordance.parquet",
        hub_42=f"{MANUSCRIPT_OUT_DIR}/3c_iatlas_instaprism_milopy_benchmark/consensus_metastate_all42_hub_summary.parquet",
        eff_42=f"{MANUSCRIPT_OUT_DIR}/3c_iatlas_instaprism_milopy_benchmark/milopy_iatlas_state_effects.parquet",
        degs_42=f"{MANUSCRIPT_OUT_DIR}/3c_iatlas_instaprism_milopy_benchmark/consensus_metastate_all42_top_degs.parquet",
    output:
        centroids=f"{COULTON_CMP_OUT_DIR}/coulton_substate_centroids.parquet",
        markers=f"{COULTON_CMP_OUT_DIR}/coulton_top_markers.parquet",
        similarity=f"{COULTON_CMP_OUT_DIR}/meta42_vs_coulton_full_similarity.parquet",
        nearest=f"{COULTON_CMP_OUT_DIR}/meta42_macrophage_nearest_coulton_substates.parquet",
        sensitivity=f"{COULTON_CMP_OUT_DIR}/meta42_coulton_sensitivity_summary.parquet",
        ici_outcomes=f"{COULTON_CMP_OUT_DIR}/coulton_substate_ici_outcomes.parquet",
        indication=f"{COULTON_CMP_OUT_DIR}/meta42_macrophage_indication_effects.parquet",
        report=f"{COULTON_CMP_OUT_DIR}/meta42_coulton_tam_comparison_report.md",
    params:
        manuscript_out_dir=MANUSCRIPT_OUT_DIR,
        out_dir=COULTON_CMP_OUT_DIR,
        batch_size=25000,
    threads: 4
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --manuscript-out-dir {params.manuscript_out_dir} \
            --output-dir {params.out_dir} \
            --batch-size {params.batch_size} \
            --skip-plots
        """


rule render_coulton_tam_comparison_figures:
    input:
        script="scripts/coulton_tam_comparison/render_figures.py",
        plotting_mod="scripts/coulton_tam_comparison/plotting.py",
        similarity=f"{COULTON_CMP_OUT_DIR}/meta42_vs_coulton_full_similarity.parquet",
        nearest=f"{COULTON_CMP_OUT_DIR}/meta42_macrophage_nearest_coulton_substates.parquet",
        ici_outcomes=f"{COULTON_CMP_OUT_DIR}/coulton_substate_ici_outcomes.parquet",
    output:
        fig1_svg=f"{COULTON_CMP_OUT_DIR}/figures/fig1_meta42_myeloid_coulton_similarity_and_ici.svg",
        fig1_png=f"{COULTON_CMP_OUT_DIR}/figures/fig1_meta42_myeloid_coulton_similarity_and_ici.png",
        fig1b_svg=f"{COULTON_CMP_OUT_DIR}/figures/fig1b_meta42_myeloid_coulton_simplified_bicluster.svg",
        fig1b_png=f"{COULTON_CMP_OUT_DIR}/figures/fig1b_meta42_myeloid_coulton_simplified_bicluster.png",
        fig2_svg=f"{COULTON_CMP_OUT_DIR}/figures/fig2_macrophage_nearest_substate_ici_forest.svg",
        fig2_png=f"{COULTON_CMP_OUT_DIR}/figures/fig2_macrophage_nearest_substate_ici_forest.png",
        fig3_svg=f"{COULTON_CMP_OUT_DIR}/figures/fig3_all42_pan_lineage_coulton_specificity.svg",
        fig3_png=f"{COULTON_CMP_OUT_DIR}/figures/fig3_all42_pan_lineage_coulton_specificity.png",
        main_cap=f"{COULTON_CMP_OUT_DIR}/coulton_tam_simplified_bicluster_main_caption.typ",
    params:
        out_dir=COULTON_CMP_OUT_DIR,
    threads: 2
    shell:
        """
        set -euo pipefail
        python {input.script} --results-dir {params.out_dir}
        """
