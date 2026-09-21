# workflow/rules/calibration.smk
# Snakemake workflow for Single-Cell Ground-Truth Calibration & Triangulation

CALIB_OUT_DIR = config.get("calib_out_dir", "output/concordance_calibration")
CALIB_RES_DIR = config.get("calib_results_dir", "results/concordance_calibration")
SC_ADATA_PATH = config.get("sc_adata", "")

CALIB_RESOLUTIONS = config.get(
    "calib_resolutions",
    "0.2 0.35 0.5 0.65 0.75 0.9 1.0 1.15 1.25 1.4 1.5 1.65 1.75 1.9 2.0 2.25 2.5 2.75 3.0 3.25 3.5"
)
CALIB_TIMEPOINTS = config.get("calib_timepoints", "Combined Pre Post")
NOISE_LEVEL = config.get("calib_noise_level", 0.05)

rule all_calibration:
    input:
        f"{CALIB_RES_DIR}/fig1_triangulation_scatter.svg",
        f"{CALIB_RES_DIR}/fig3_resolution_concordance_curve.svg",
        f"{CALIB_RES_DIR}/fig4_discrepancy_attribution_bar.svg",
        f"{CALIB_OUT_DIR}/triangulation_diagnostics_summary.parquet",
        f"{CALIB_OUT_DIR}/resolution_trend_summary.parquet"

rule sc_extract_ground_truth:
    input:
        script="scripts/concordance_calibration/01_extract_sc_ground_truth.py"
    output:
        counts=f"{CALIB_OUT_DIR}/sc_ground_truth_counts.parquet",
        umi=f"{CALIB_OUT_DIR}/sc_ground_truth_umi.parquet"
    params:
        adata=SC_ADATA_PATH,
        timepoints=CALIB_TIMEPOINTS,
        resolutions=CALIB_RESOLUTIONS,
        out_dir=CALIB_OUT_DIR
    log:
        f"{CALIB_OUT_DIR}/logs/01_extract_ground_truth.log"
    shell:
        """
        python {input.script} \
            --timepoints {params.timepoints} \
            --resolutions {params.resolutions} \
            --out-dir {params.out_dir} > {log} 2>&1
        """

rule sc_pseudobulk_self_deconvolution:
    input:
        script="scripts/concordance_calibration/02_pseudobulk_self_deconvolution.py",
        counts=f"{CALIB_OUT_DIR}/sc_ground_truth_counts.parquet",
        umi=f"{CALIB_OUT_DIR}/sc_ground_truth_umi.parquet"
    output:
        deconv=f"{CALIB_OUT_DIR}/sc_pseudobulk_deconv_fractions.parquet",
        fidelity=f"{CALIB_OUT_DIR}/deconv_fidelity_benchmark.parquet"
    params:
        noise_level=NOISE_LEVEL,
        out_dir=CALIB_OUT_DIR
    log:
        f"{CALIB_OUT_DIR}/logs/02_pseudobulk_deconv.log"
    shell:
        """
        python {input.script} \
            --gt-counts {input.counts} \
            --gt-umi {input.umi} \
            --noise-level {params.noise_level} \
            --out-dir {params.out_dir} > {log} 2>&1
        """

rule sc_fit_triangulation_models:
    input:
        script="scripts/concordance_calibration/03_fit_triangulation_models.py",
        counts=f"{CALIB_OUT_DIR}/sc_ground_truth_counts.parquet",
        umi=f"{CALIB_OUT_DIR}/sc_ground_truth_umi.parquet",
        deconv=f"{CALIB_OUT_DIR}/sc_pseudobulk_deconv_fractions.parquet"
    output:
        triangulation=f"{CALIB_OUT_DIR}/triangulation_effect_estimates.parquet"
    params:
        out_dir=CALIB_OUT_DIR
    log:
        f"{CALIB_OUT_DIR}/logs/03_fit_triangulation.log"
    shell:
        """
        python {input.script} \
            --counts {input.counts} \
            --umi {input.umi} \
            --deconv {input.deconv} \
            --out-dir {params.out_dir} > {log} 2>&1
        """

rule sc_diagnose_discrepancies:
    input:
        script="scripts/concordance_calibration/04_diagnose_discrepancies.py",
        triangulation=f"{CALIB_OUT_DIR}/triangulation_effect_estimates.parquet",
        fidelity=f"{CALIB_OUT_DIR}/deconv_fidelity_benchmark.parquet"
    output:
        diagnostics=f"{CALIB_OUT_DIR}/triangulation_diagnostics_summary.parquet",
        trends=f"{CALIB_OUT_DIR}/resolution_trend_summary.parquet"
    params:
        out_dir=CALIB_OUT_DIR
    log:
        f"{CALIB_OUT_DIR}/logs/04_diagnose_discrepancies.log"
    shell:
        """
        python {input.script} \
            --triangulation {input.triangulation} \
            --fidelity {input.fidelity} \
            --out-dir {params.out_dir} > {log} 2>&1
        """

rule sc_plot_triangulation_figures:
    input:
        script="scripts/concordance_calibration/05_plot_triangulation_figures.py",
        diagnostics=f"{CALIB_OUT_DIR}/triangulation_diagnostics_summary.parquet",
        trends=f"{CALIB_OUT_DIR}/resolution_trend_summary.parquet",
        fidelity=f"{CALIB_OUT_DIR}/deconv_fidelity_benchmark.parquet"
    output:
        fig1=f"{CALIB_RES_DIR}/fig1_triangulation_scatter.svg",
        fig3=f"{CALIB_RES_DIR}/fig3_resolution_concordance_curve.svg",
        fig4=f"{CALIB_RES_DIR}/fig4_discrepancy_attribution_bar.svg"
    params:
        results_dir=CALIB_RES_DIR
    log:
        f"{CALIB_RES_DIR}/logs/05_plot_figures.log"
    shell:
        """
        python {input.script} \
            --diagnostics {input.diagnostics} \
            --trends {input.trends} \
            --fidelity {input.fidelity} \
            --results-dir {params.results_dir} > {log} 2>&1
        """
