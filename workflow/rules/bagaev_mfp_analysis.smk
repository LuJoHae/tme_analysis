"""Snakemake workflow for Bagaev et al. (2021) MFP TME Subtyping, Deconvolution & Single-Cell Pseudobulk Analysis."""

import os

DATA_DIR = os.environ.get("TME_DATA_DIR", "data")
BAGAEV_OUT_DIR = f"{DATA_DIR}/output/bagaev_mfp_analysis"
BAGAEV_RESULTS_DIR = f"{DATA_DIR}/results/bagaev_mfp_analysis"
LAIR_DIR = "/storage/halu/lair" if os.path.exists("/storage/halu/lair") else f"{DATA_DIR}/lair"
BENCH_3C_DIR = (
    f"{DATA_DIR}/output/manuscript_pipeline/3c_iatlas_instaprism_milopy_benchmark"
    if os.path.exists(f"{DATA_DIR}/output/manuscript_pipeline/3c_iatlas_instaprism_milopy_benchmark/oof_9a_9b_patient_landscape_and_curves.parquet")
    else "output/manuscript_pipeline/3c_iatlas_instaprism_milopy_benchmark"
)
SC9_DECONV_DIR = (
    f"{DATA_DIR}/output/sc9_milo_deconv_validation/02_deconvolution"
    if os.path.exists(f"{DATA_DIR}/output/sc9_milo_deconv_validation/02_deconvolution")
    else "output/sc9_milo_deconv_validation/02_deconvolution"
)
SC9_MILO_DIR = (
    f"{DATA_DIR}/output/sc9_milo_deconv_validation/04_milopy"
    if os.path.exists(f"{DATA_DIR}/output/sc9_milo_deconv_validation/04_milopy")
    else "output/sc9_milo_deconv_validation/04_milopy"
)


rule all_bagaev_mfp_analysis:
    input:
        f"{BAGAEV_OUT_DIR}/iatlas_bagaev_sample_subtypes.parquet",
        f"{BAGAEV_OUT_DIR}/iatlas_subtype_response_rates.parquet",
        f"{BAGAEV_OUT_DIR}/iatlas_meta42_enrichment_by_bagaev_subtype.parquet",
        f"{BAGAEV_OUT_DIR}/iatlas_within_subtype_deconv_auc.parquet",
        f"{BAGAEV_OUT_DIR}/sc23_pseudobulk_cohort_summary.parquet",
        f"{BAGAEV_OUT_DIR}/sc9_milopy_subtype_response_rates.parquet",
        f"{BAGAEV_RESULTS_DIR}/fig1_iatlas_29fges_portrait_and_clustering.svg",
        f"{BAGAEV_RESULTS_DIR}/fig2_iatlas_subtype_response_rates.svg",
        f"{BAGAEV_RESULTS_DIR}/fig3_iatlas_deconvolution_integration.svg",
        f"{BAGAEV_RESULTS_DIR}/fig4_sc23_pseudobulk_clustering_conservation.svg",


rule bagaev_01_classify_and_cluster_iatlas:
    input:
        script="scripts/bagaev_mfp_analysis/01_classify_and_cluster_iatlas.py",
    output:
        subtypes=f"{BAGAEV_OUT_DIR}/iatlas_bagaev_sample_subtypes.parquet",
        resp=f"{BAGAEV_OUT_DIR}/iatlas_subtype_response_rates.parquet",
        tcga_ref=f"{BAGAEV_OUT_DIR}/tcga_bagaev_reference_29fges.parquet",
    params:
        out_dir=BAGAEV_OUT_DIR,
        lair_dir=LAIR_DIR,
    threads: 4
    shell:
        "python3 {input.script} --out-dir {params.out_dir} --lair-dir {params.lair_dir}"


rule bagaev_02_integrate_iatlas_deconvolution:
    input:
        script="scripts/bagaev_mfp_analysis/02_integrate_iatlas_deconvolution.py",
        subtypes=f"{BAGAEV_OUT_DIR}/iatlas_bagaev_sample_subtypes.parquet",
    output:
        enr=f"{BAGAEV_OUT_DIR}/iatlas_meta42_enrichment_by_bagaev_subtype.parquet",
        within=f"{BAGAEV_OUT_DIR}/iatlas_within_subtype_deconv_auc.parquet",
        cv=f"{BAGAEV_OUT_DIR}/iatlas_bagaev_vs_deconv_cv_benchmark.parquet",
    params:
        bagaev_dir=BAGAEV_OUT_DIR,
        bench_3c_dir=BENCH_3C_DIR,
        sc9_deconv_dir=SC9_DECONV_DIR,
    threads: 4
    shell:
        "python3 {input.script} --bagaev-dir {params.bagaev_dir} --bench-3c-dir {params.bench_3c_dir} --sc9-deconv-dir {params.sc9_deconv_dir}"


rule bagaev_03_analyze_sc_pseudobulk_and_milopy:
    input:
        script="scripts/bagaev_mfp_analysis/03_analyze_sc_pseudobulk_and_milopy.py",
        tcga_ref=f"{BAGAEV_OUT_DIR}/tcga_bagaev_reference_29fges.parquet",
    output:
        summary=f"{BAGAEV_OUT_DIR}/sc23_pseudobulk_cohort_summary.parquet",
        samples=f"{BAGAEV_OUT_DIR}/sc23_pseudobulk_bagaev_samples.parquet",
        sc9_resp=f"{BAGAEV_OUT_DIR}/sc9_milopy_subtype_response_rates.parquet",
    params:
        out_dir=BAGAEV_OUT_DIR,
        preprocessed_dir=f"{DATA_DIR}/preprocessed",
        lair_dir=LAIR_DIR,
        pb_cache_dir=f"{DATA_DIR}/pseudobulk",
        sc9_milo_dir=SC9_MILO_DIR,
    threads: 8
    shell:
        "python3 {input.script} --out-dir {params.out_dir} --preprocessed-dir {params.preprocessed_dir} --lair-dir {params.lair_dir} --pb-cache-dir {params.pb_cache_dir} --sc9-milo-dir {params.sc9_milo_dir}"


rule bagaev_04_plot_figures:
    input:
        script="scripts/bagaev_mfp_analysis/04_plot_bagaev_mfp_figures.py",
        subtypes=f"{BAGAEV_OUT_DIR}/iatlas_bagaev_sample_subtypes.parquet",
        resp=f"{BAGAEV_OUT_DIR}/iatlas_subtype_response_rates.parquet",
        enr=f"{BAGAEV_OUT_DIR}/iatlas_meta42_enrichment_by_bagaev_subtype.parquet",
        within=f"{BAGAEV_OUT_DIR}/iatlas_within_subtype_deconv_auc.parquet",
        sc_sum=f"{BAGAEV_OUT_DIR}/sc23_pseudobulk_cohort_summary.parquet",
        sc9_resp=f"{BAGAEV_OUT_DIR}/sc9_milopy_subtype_response_rates.parquet",
    output:
        fig1=f"{BAGAEV_RESULTS_DIR}/fig1_iatlas_29fges_portrait_and_clustering.svg",
        fig2=f"{BAGAEV_RESULTS_DIR}/fig2_iatlas_subtype_response_rates.svg",
        fig3=f"{BAGAEV_RESULTS_DIR}/fig3_iatlas_deconvolution_integration.svg",
        fig4=f"{BAGAEV_RESULTS_DIR}/fig4_sc23_pseudobulk_clustering_conservation.svg",
    params:
        data_dir=BAGAEV_OUT_DIR,
        results_dir=BAGAEV_RESULTS_DIR,
    threads: 2
    shell:
        "python3 {input.script} --data-dir {params.data_dir} --results-dir {params.results_dir}"
