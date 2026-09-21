# workflow/rules/concordance.smk
# Snakemake workflow for Single-Cell and Bulk Deconvolution Concordance Analysis

CONCORDANCE_OUT_DIR = config.get("concordance_out_dir", "output/concordance")
CONCORDANCE_RES_DIR = config.get("concordance_results_dir", "results/concordance")
SC_ADATA_PATH = config.get("sc_adata", "/storage/halu/data/GSE120575/gse120575_processed.h5ad")
BULK_COUNTS_PATH = config.get("bulk_counts", "data/bulk_counts.parquet")
CLINICAL_PATH = config.get("clinical_meta", "data/bulk_clinical.parquet")

PATIENT_COL = config.get("sc_patient_col", "patient")
RESPONSE_COL = config.get("sc_response_col", "response")
CLUSTER_COL = config.get("sc_cluster_col", "celltypist_leiden_0.5")
LINEAGE_COL = config.get("sc_lineage_col", "cell_type")
MALIGNANT_LABEL = config.get("malignant_label", "Malignant")
BULK_SAMPLE_COL = config.get("bulk_sample_col", "sample_id")
BULK_RESPONSE_COL = config.get("bulk_response_col", "response")

rule all_concordance:
    input:
        f"{CONCORDANCE_RES_DIR}/fig1_concordance_scatter.svg",
        f"{CONCORDANCE_RES_DIR}/fig2_purity_adjustment_shift.svg",
        f"{CONCORDANCE_RES_DIR}/fig3_forest_effect_sizes.svg",
        f"{CONCORDANCE_RES_DIR}/fig4_insilico_calibration.svg",
        f"{CONCORDANCE_RES_DIR}/fig6_hierarchical_resolution_comparison.svg",
        f"{CONCORDANCE_OUT_DIR}/concordance_metrics_summary.parquet",
        f"{CONCORDANCE_OUT_DIR}/concordance_diagnostics.parquet"

rule sc_patient_response_association:
    input:
        script="scripts/concordance/01_sc_patient_response_association.py",
        adata=SC_ADATA_PATH
    output:
        effects=f"{CONCORDANCE_OUT_DIR}/sc_patient_response_effects.parquet"
    params:
        patient_col=PATIENT_COL,
        response_col=RESPONSE_COL,
        cluster_col=CLUSTER_COL
    log:
        f"{CONCORDANCE_OUT_DIR}/logs/01_sc_association.log"
    shell:
        """
        python {input.script} \
            --adata {input.adata} \
            --patient-col {params.patient_col} \
            --response-col {params.response_col} \
            --cluster-col {params.cluster_col} \
            --out-parquet {output.effects} > {log} 2>&1
        """

rule build_curated_reference:
    input:
        script="scripts/concordance/02_build_curated_reference.py",
        adata=SC_ADATA_PATH
    output:
        phi=f"{CONCORDANCE_OUT_DIR}/curated_reference_phi.parquet",
        mrna=f"{CONCORDANCE_OUT_DIR}/mrna_scaling_factors.parquet",
        collinearity=f"{CONCORDANCE_OUT_DIR}/reference_collinearity.parquet",
        hierarchy=f"{CONCORDANCE_OUT_DIR}/cell_hierarchy.parquet"
    params:
        cluster_col=CLUSTER_COL,
        malignant_label=MALIGNANT_LABEL,
        out_dir=CONCORDANCE_OUT_DIR
    log:
        f"{CONCORDANCE_OUT_DIR}/logs/02_build_reference.log"
    shell:
        """
        python {input.script} \
            --adata {input.adata} \
            --cluster-col {params.cluster_col} \
            --malignant-label {params.malignant_label} \
            --out-dir {params.out_dir} > {log} 2>&1
        """

rule insilico_pseudobulk_validation:
    input:
        script="scripts/concordance/03_insilico_pseudobulk_validation.py",
        adata=SC_ADATA_PATH,
        reference=f"{CONCORDANCE_OUT_DIR}/curated_reference_phi.parquet",
        hierarchy=f"{CONCORDANCE_OUT_DIR}/cell_hierarchy.parquet",
        sc_effects=f"{CONCORDANCE_OUT_DIR}/sc_patient_response_effects.parquet"
    output:
        metrics=f"{CONCORDANCE_OUT_DIR}/insilico_validation_metrics.parquet",
        comparison=f"{CONCORDANCE_OUT_DIR}/insilico_recovery_comparison.parquet",
        fractions=f"{CONCORDANCE_OUT_DIR}/insilico_pseudobulk_fractions.parquet"
    params:
        patient_col=PATIENT_COL,
        cluster_col=CLUSTER_COL,
        response_col=RESPONSE_COL,
        out_dir=CONCORDANCE_OUT_DIR
    log:
        f"{CONCORDANCE_OUT_DIR}/logs/03_insilico_validation.log"
    shell:
        """
        python {input.script} \
            --adata {input.adata} \
            --patient-col {params.patient_col} \
            --cluster-col {params.cluster_col} \
            --response-col {params.response_col} \
            --reference {input.reference} \
            --hierarchy {input.hierarchy} \
            --sc-effects {input.sc_effects} \
            --out-dir {params.out_dir} > {log} 2>&1
        """

rule deconvolve_bulk_bayesprism:
    input:
        script="scripts/concordance/04_deconvolve_bulk_bayesprism.py",
        bulk=BULK_COUNTS_PATH,
        reference=f"{CONCORDANCE_OUT_DIR}/curated_reference_phi.parquet",
        hierarchy=f"{CONCORDANCE_OUT_DIR}/cell_hierarchy.parquet",
        mrna=f"{CONCORDANCE_OUT_DIR}/mrna_scaling_factors.parquet"
    output:
        raw=f"{CONCORDANCE_OUT_DIR}/bulk_deconv_fractions_raw.parquet",
        norm=f"{CONCORDANCE_OUT_DIR}/bulk_deconv_fractions_normalized.parquet",
        mrna_scaled=f"{CONCORDANCE_OUT_DIR}/bulk_deconv_fractions_mrna_scaled.parquet"
    params:
        malignant_label=MALIGNANT_LABEL,
        out_dir=CONCORDANCE_OUT_DIR
    threads: 4
    log:
        f"{CONCORDANCE_OUT_DIR}/logs/04_deconvolve_bulk.log"
    shell:
        """
        python {input.script} \
            --bulk-counts {input.bulk} \
            --reference {input.reference} \
            --hierarchy {input.hierarchy} \
            --mrna-scaling {input.mrna} \
            --malignant-label {params.malignant_label} \
            --out-dir {params.out_dir} > {log} 2>&1
        """

rule bulk_response_association:
    input:
        script="scripts/concordance/05_bulk_response_association.py",
        raw=f"{CONCORDANCE_OUT_DIR}/bulk_deconv_fractions_raw.parquet",
        norm=f"{CONCORDANCE_OUT_DIR}/bulk_deconv_fractions_normalized.parquet",
        mrna=f"{CONCORDANCE_OUT_DIR}/bulk_deconv_fractions_mrna_scaled.parquet",
        clinical=CLINICAL_PATH
    output:
        effects=f"{CONCORDANCE_OUT_DIR}/bulk_response_effects.parquet"
    params:
        sample_col=BULK_SAMPLE_COL,
        response_col=BULK_RESPONSE_COL
    log:
        f"{CONCORDANCE_OUT_DIR}/logs/05_bulk_association.log"
    shell:
        """
        python {input.script} \
            --fractions-raw {input.raw} \
            --fractions-norm {input.norm} \
            --fractions-mrna {input.mrna} \
            --clinical {input.clinical} \
            --sample-col {params.sample_col} \
            --response-col {params.response_col} \
            --out-parquet {output.effects} > {log} 2>&1
        """

rule evaluate_concordance:
    input:
        script="scripts/concordance/06_evaluate_concordance.py",
        sc_effects=f"{CONCORDANCE_OUT_DIR}/sc_patient_response_effects.parquet",
        bulk_effects=f"{CONCORDANCE_OUT_DIR}/bulk_response_effects.parquet",
        insilico=f"{CONCORDANCE_OUT_DIR}/insilico_recovery_comparison.parquet",
        collinearity=f"{CONCORDANCE_OUT_DIR}/reference_collinearity.parquet"
    output:
        summary=f"{CONCORDANCE_OUT_DIR}/concordance_metrics_summary.parquet",
        details=f"{CONCORDANCE_OUT_DIR}/concordance_state_details.parquet",
        diagnostics=f"{CONCORDANCE_OUT_DIR}/concordance_diagnostics.parquet"
    params:
        out_dir=CONCORDANCE_OUT_DIR
    log:
        f"{CONCORDANCE_OUT_DIR}/logs/06_evaluate_concordance.log"
    shell:
        """
        python {input.script} \
            --sc-effects {input.sc_effects} \
            --bulk-effects {input.bulk_effects} \
            --insilico {input.insilico} \
            --collinearity {input.collinearity} \
            --out-dir {params.out_dir} > {log} 2>&1
        """

rule plot_concordance_figures:
    input:
        script="scripts/concordance/07_plot_concordance_figures.py",
        summary=f"{CONCORDANCE_OUT_DIR}/concordance_metrics_summary.parquet",
        details=f"{CONCORDANCE_OUT_DIR}/concordance_state_details.parquet",
        insilico=f"{CONCORDANCE_OUT_DIR}/insilico_recovery_comparison.parquet"
    output:
        fig1=f"{CONCORDANCE_RES_DIR}/fig1_concordance_scatter.svg",
        fig2=f"{CONCORDANCE_RES_DIR}/fig2_purity_adjustment_shift.svg",
        fig3=f"{CONCORDANCE_RES_DIR}/fig3_forest_effect_sizes.svg",
        fig4=f"{CONCORDANCE_RES_DIR}/fig4_insilico_calibration.svg"
    params:
        results_dir=CONCORDANCE_RES_DIR
    log:
        f"{CONCORDANCE_RES_DIR}/logs/07_plot_figures.log"
    shell:
        """
        python {input.script} \
            --metrics-summary {input.summary} \
            --state-details {input.details} \
            --insilico {input.insilico} \
            --results-dir {params.results_dir} > {log} 2>&1
        """

rule hierarchical_bayesprism_resolution:
    input:
        script="scripts/concordance/09_hierarchical_bayesprism_resolution.py",
        ref="output/output/sade_feldman_deconv_validation/reference_phi_res0.5.parquet",
        sc_effects=f"{CONCORDANCE_OUT_DIR}/sc_patient_response_effects.parquet"
    output:
        effects=f"{CONCORDANCE_OUT_DIR}/hierarchical_bayesprism_effects.parquet",
        metrics=f"{CONCORDANCE_OUT_DIR}/hierarchical_concordance_metrics.parquet",
        fig6=f"{CONCORDANCE_RES_DIR}/fig6_hierarchical_resolution_comparison.svg"
    params:
        out_dir=CONCORDANCE_OUT_DIR,
        results_dir=CONCORDANCE_RES_DIR,
        threshold=0.85
    log:
        f"{CONCORDANCE_RES_DIR}/logs/09_hierarchical_resolution.log"
    shell:
        """
        python {input.script} \
            --ref {input.ref} \
            --sc-effects {input.sc_effects} \
            --threshold {params.threshold} \
            --out-dir {params.out_dir} \
            --results-dir {params.results_dir} > {log} 2>&1
        """
