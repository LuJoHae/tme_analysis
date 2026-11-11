# workflow/rules/sc9_deconv_validation.smk
# Snakemake workflow module for the 9-Cohort Multi-Resolution Single-Cell Reference
# Deconvolution & Milo Validation Pipeline.
# Compatible with both OLM (/storage/halu) and CSCS Alps Eiger (/ritom/scratch/cscs/lhaeuser/data-test).

import os

SC9_PIPELINE_ROOT = os.environ.get(
    "SC9_PIPELINE_OUTPUT_ROOT",
    os.path.join(DATA_DIR, "output/sc9_milo_deconv_validation"),
)
SC9_RESULTS_DIR = os.environ.get(
    "SC9_RESULTS_DIR",
    os.path.join(DATA_DIR, "results/sc9_milo_deconv_validation"),
)

SC9_DIR_REF = f"{SC9_PIPELINE_ROOT}/01_references"
SC9_DIR_DECONV = f"{SC9_PIPELINE_ROOT}/02_deconvolution"
SC9_DIR_LR = f"{SC9_PIPELINE_ROOT}/03_logistic_regression"
SC9_DIR_MILO = f"{SC9_PIPELINE_ROOT}/04_milopy"
SC9_DIR_CONC = f"{SC9_PIPELINE_ROOT}/05_concordance"

SC9_ALL_COHORTS = [
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

SC9_RESOLUTIONS = ["0.25", "0.5", "0.75", "1.0", "1.25", "1.5", "1.75", "2.0"]


rule sc9_build_single_reference:
    input:
        script="scripts/sc9_milo_deconv_validation/01_build_single_reference.py",
        adata=f"{DATA_DIR}/preprocessed/{{sc9_cohort}}.h5ad",
    output:
        phi_05=f"{SC9_DIR_REF}/reference_phi_{{sc9_cohort}}_res0.5.parquet",
        metrics=f"{SC9_DIR_REF}/reference_resolution_metrics_{{sc9_cohort}}.parquet",
    wildcard_constraints:
        sc9_cohort="GSE120575|CELLxGENE_7b20c613|CELLxGENE_05a8c945|CELLxGENE_6f9de485|GSE207422|GSE243013|GSE233203|GSE200996|GSE316195"
    params:
        out_dir=SC9_DIR_REF,
    shell:
        """
        python {input.script} --cohort {wildcards.sc9_cohort} --adata {input.adata} --out-dir {params.out_dir}
        """


rule sc9_build_combined_references:
    input:
        script="scripts/sc9_milo_deconv_validation/02_build_combined_references.py",
        single_metrics=expand(
            f"{SC9_DIR_REF}/reference_resolution_metrics_{{sc9_cohort}}.parquet",
            sc9_cohort=SC9_ALL_COHORTS,
        ),
    output:
        all9_phi_05=f"{SC9_DIR_REF}/reference_phi_all_9_combined_res0.5.parquet",
        all9_metrics=f"{SC9_DIR_REF}/reference_resolution_metrics_all_9_combined.parquet",
    params:
        out_dir=SC9_DIR_REF,
    shell:
        """
        python {input.script} --combo all --out-dir {params.out_dir}
        """


rule sc9_deconvolute_iatlas:
    input:
        script="scripts/sc9_milo_deconv_validation/03_deconvolute_iatlas.py",
        all9_metrics=f"{SC9_DIR_REF}/reference_resolution_metrics_all_9_combined.parquet",
        single_metrics=expand(
            f"{SC9_DIR_REF}/reference_resolution_metrics_{{sc9_cohort}}.parquet",
            sc9_cohort=SC9_ALL_COHORTS,
        ),
    output:
        touch_done=f"{SC9_DIR_DECONV}/deconv_completed.flag",
        fracs_all9_05=f"{SC9_DIR_DECONV}/deconv_fractions_all_9_combined_res0.5.parquet",
    params:
        ref_dir=SC9_DIR_REF,
        out_dir=SC9_DIR_DECONV,
    shell:
        """
        python {input.script} --ref-dir {params.ref_dir} --out-dir {params.out_dir}
        touch {output.touch_done}
        """


rule sc9_logistic_regression:
    input:
        script="scripts/sc9_milo_deconv_validation/04_logistic_regression.py",
        touch_done=f"{SC9_DIR_DECONV}/deconv_completed.flag",
        fracs_all9_05=f"{SC9_DIR_DECONV}/deconv_fractions_all_9_combined_res0.5.parquet",
    output:
        results=f"{SC9_DIR_LR}/logistic_regression_results_all_resolutions.parquet",
    params:
        deconv_dir=SC9_DIR_DECONV,
        out_dir=SC9_DIR_LR,
    shell:
        """
        python {input.script} --deconv-dir {params.deconv_dir} --out-dir {params.out_dir}
        """


rule sc9_analyze_milopy:
    input:
        script="scripts/sc9_milo_deconv_validation/05_analyze_milopy.py",
        adata=f"{DATA_DIR}/preprocessed/GSE120575.h5ad",
    output:
        milo_results=f"{SC9_DIR_MILO}/milopy_cell_state_da_all_resolutions.parquet",
    params:
        out_dir=SC9_DIR_MILO,
    shell:
        """
        python {input.script} --cohort GSE120575 --out-dir {params.out_dir}
        """


rule sc9_compare_concordance:
    input:
        script="scripts/sc9_milo_deconv_validation/06_compare_concordance.py",
        logistic=f"{SC9_DIR_LR}/logistic_regression_results_all_resolutions.parquet",
        milo=f"{SC9_DIR_MILO}/milopy_cell_state_da_all_resolutions.parquet",
        all9_metrics=f"{SC9_DIR_REF}/reference_resolution_metrics_all_9_combined.parquet",
    output:
        bench_summary=f"{SC9_DIR_CONC}/multi_resolution_benchmark_summary.parquet",
    params:
        ref_dir=SC9_DIR_REF,
        logistic_dir=SC9_DIR_LR,
        milo_dir=SC9_DIR_MILO,
        out_dir=SC9_DIR_CONC,
    shell:
        """
        python {input.script} --ref-dir {params.ref_dir} --logistic-dir {params.logistic_dir} --milo-dir {params.milo_dir} --out-dir {params.out_dir}
        """


rule sc9_plot_validation_figures:
    input:
        script="scripts/sc9_milo_deconv_validation/07_plot_figures.py",
        bench_summary=f"{SC9_DIR_CONC}/multi_resolution_benchmark_summary.parquet",
    output:
        fig1_svg=f"{SC9_RESULTS_DIR}/fig1_predictive_capacity_benchmark.svg",
        fig1_png=f"{SC9_RESULTS_DIR}/fig1_predictive_capacity_benchmark.png",
        fig2_svg=f"{SC9_RESULTS_DIR}/fig2_signature_collinearity_benchmark.svg",
        fig3_svg=f"{SC9_RESULTS_DIR}/fig3_reference_selection_decision_framework.svg",
        fig4_svg=f"{SC9_RESULTS_DIR}/fig4_winning_classifier_biomarker_drivers.svg",
        fig5_svg=f"{SC9_RESULTS_DIR}/fig5_selected_cohort_auc_and_marginal_gain.svg",
    params:
        data_dir=SC9_DIR_CONC,
        results_dir=SC9_RESULTS_DIR,
    shell:
        """
        python {input.script} --data-dir {params.data_dir} --results-dir {params.results_dir}
        """


rule run_sc9_deconv_pipeline:
    input:
        f"{SC9_RESULTS_DIR}/fig1_predictive_capacity_benchmark.svg",
        f"{SC9_RESULTS_DIR}/fig1_predictive_capacity_benchmark.png",
        f"{SC9_RESULTS_DIR}/fig2_signature_collinearity_benchmark.svg",
        f"{SC9_RESULTS_DIR}/fig3_reference_selection_decision_framework.svg",
        f"{SC9_RESULTS_DIR}/fig4_winning_classifier_biomarker_drivers.svg",
        f"{SC9_RESULTS_DIR}/fig5_selected_cohort_auc_and_marginal_gain.svg",
