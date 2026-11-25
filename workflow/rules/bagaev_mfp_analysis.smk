"""Snakemake workflow for Bagaev et al. (2021) MFP TME Subtyping, Deconvolution & Single-Cell Pseudobulk Analysis."""

import os

DATA_DIR = os.environ.get("TME_DATA_DIR", "data")
LAIR_DIR = os.environ.get("TME_LAIR_DIR", f"{DATA_DIR}/lair")
BAGAEV_OUT_DIR = f"{DATA_DIR}/output/bagaev_mfp_analysis"
BAGAEV_RESULTS_DIR = f"{DATA_DIR}/results/bagaev_mfp_analysis"
BENCH_3C_DIR = f"{DATA_DIR}/output/manuscript_pipeline/3c_iatlas_instaprism_milopy_benchmark"
SC9_DECONV_DIR = f"{DATA_DIR}/output/sc9_milo_deconv_validation/02_deconvolution"
SC9_MILO_DIR = f"{DATA_DIR}/output/sc9_milo_deconv_validation/04_milopy"


rule all_bagaev_mfp_analysis:
    input:
        f"{BAGAEV_OUT_DIR}/iatlas_bagaev_sample_subtypes.parquet",
        f"{BAGAEV_OUT_DIR}/iatlas_subtype_29fges_means.parquet",
        f"{BAGAEV_OUT_DIR}/iatlas_subtype_response_rates.parquet",
        f"{BAGAEV_OUT_DIR}/iatlas_meta42_enrichment_by_bagaev_subtype.parquet",
        f"{BAGAEV_OUT_DIR}/iatlas_meta42_enrichment_heatmap_z.parquet",
        f"{BAGAEV_OUT_DIR}/iatlas_deconv_subtype_correlations.parquet",
        f"{BAGAEV_OUT_DIR}/iatlas_deconv_cluster_response_rates.parquet",
        f"{BAGAEV_OUT_DIR}/iatlas_multires_deconv_cluster_benchmark.parquet",
        f"{BAGAEV_OUT_DIR}/iatlas_within_subtype_deconv_auc.parquet",
        f"{BAGAEV_OUT_DIR}/iatlas_bagaev_vs_deconv_cv_benchmark.parquet",
        f"{BAGAEV_OUT_DIR}/iatlas_bagaev_vs_deconv_per_cohort_cv.parquet",
        f"{BAGAEV_OUT_DIR}/sc23_pseudobulk_cohort_summary.parquet",
        f"{BAGAEV_OUT_DIR}/sc9_milopy_subtype_response_rates.parquet",
        f"{BAGAEV_OUT_DIR}/sc9_milopy_cohort_bagaev_aucs.parquet",
        f"{BAGAEV_RESULTS_DIR}/fig1_iatlas_29fges_portrait_and_clustering.svg",
        f"{BAGAEV_RESULTS_DIR}/fig2_iatlas_subtype_response_rates.svg",
        f"{BAGAEV_RESULTS_DIR}/fig3_iatlas_deconvolution_integration.svg",
        f"{BAGAEV_RESULTS_DIR}/fig4_sc23_pseudobulk_clustering_conservation.svg",
        f"{BAGAEV_RESULTS_DIR}/fig5_sc9_milopy_response_correlation.svg",
        f"{BAGAEV_RESULTS_DIR}/fig6_deconv_subtype_correlations.svg",
        f"{BAGAEV_RESULTS_DIR}/fig7_deconv_sample_clusters_and_response.svg",


rule bagaev_01_classify_and_cluster_iatlas:
    input:
        script="scripts/bagaev_mfp_analysis/01_classify_and_cluster_iatlas.py",
    output:
        subtypes=f"{BAGAEV_OUT_DIR}/iatlas_bagaev_sample_subtypes.parquet",
        means=f"{BAGAEV_OUT_DIR}/iatlas_subtype_29fges_means.parquet",
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
        raw_frac=f"{BENCH_3C_DIR}/fractions/raw_instaprism_Consensus-Meta42_on_PanCancer-Combined.parquet",
        hub=f"{BENCH_3C_DIR}/consensus_metastate_all42_hub_summary.parquet",
        degs=f"{BENCH_3C_DIR}/consensus_metastate_all42_top_degs.parquet",
    output:
        enr=f"{BAGAEV_OUT_DIR}/iatlas_meta42_enrichment_by_bagaev_subtype.parquet",
        hm_z=f"{BAGAEV_OUT_DIR}/iatlas_meta42_enrichment_heatmap_z.parquet",
        joined=f"{BAGAEV_OUT_DIR}/iatlas_bagaev_joined_meta42_oof.parquet",
        sub_corr=f"{BAGAEV_OUT_DIR}/iatlas_deconv_subtype_correlations.parquet",
        fges_corr=f"{BAGAEV_OUT_DIR}/iatlas_meta42_vs_29fges_correlations.parquet",
        cl_samp=f"{BAGAEV_OUT_DIR}/iatlas_deconv_cluster_samples.parquet",
        cl_cent=f"{BAGAEV_OUT_DIR}/iatlas_deconv_cluster_centroids_k4.parquet",
        cl_resp=f"{BAGAEV_OUT_DIR}/iatlas_deconv_cluster_response_rates.parquet",
        cl_cont=f"{BAGAEV_OUT_DIR}/iatlas_deconv_clusters_vs_bagaev_contingency.parquet",
        multires=f"{BAGAEV_OUT_DIR}/iatlas_multires_deconv_cluster_benchmark.parquet",
        within=f"{BAGAEV_OUT_DIR}/iatlas_within_subtype_deconv_auc.parquet",
        cv=f"{BAGAEV_OUT_DIR}/iatlas_bagaev_vs_deconv_cv_benchmark.parquet",
        cv_coh=f"{BAGAEV_OUT_DIR}/iatlas_bagaev_vs_deconv_per_cohort_cv.parquet",
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
        sc9_auc=f"{BAGAEV_OUT_DIR}/sc9_milopy_cohort_bagaev_aucs.parquet",
    params:
        out_dir=BAGAEV_OUT_DIR,
        preprocessed_dir=f"{DATA_DIR}/preprocessed",
        ensembl_dir=f"{DATA_DIR}/ensembl",
        lair_dir=LAIR_DIR,
        pb_cache_dir=f"{DATA_DIR}/pseudobulk",
        sc9_milo_dir=SC9_MILO_DIR,
    threads: 8
    shell:
        "python3 {input.script} --out-dir {params.out_dir} --preprocessed-dir {params.preprocessed_dir} --ensembl-dir {params.ensembl_dir} --lair-dir {params.lair_dir} --pb-cache-dir {params.pb_cache_dir} --sc9-milo-dir {params.sc9_milo_dir}"


rule bagaev_04_plot_figures:
    input:
        script="scripts/bagaev_mfp_analysis/04_plot_bagaev_mfp_figures.py",
        subtypes=f"{BAGAEV_OUT_DIR}/iatlas_bagaev_sample_subtypes.parquet",
        means=f"{BAGAEV_OUT_DIR}/iatlas_subtype_29fges_means.parquet",
        resp=f"{BAGAEV_OUT_DIR}/iatlas_subtype_response_rates.parquet",
        enr=f"{BAGAEV_OUT_DIR}/iatlas_meta42_enrichment_by_bagaev_subtype.parquet",
        hm_z=f"{BAGAEV_OUT_DIR}/iatlas_meta42_enrichment_heatmap_z.parquet",
        joined=f"{BAGAEV_OUT_DIR}/iatlas_bagaev_joined_meta42_oof.parquet",
        sub_corr=f"{BAGAEV_OUT_DIR}/iatlas_deconv_subtype_correlations.parquet",
        fges_corr=f"{BAGAEV_OUT_DIR}/iatlas_meta42_vs_29fges_correlations.parquet",
        cl_samp=f"{BAGAEV_OUT_DIR}/iatlas_deconv_cluster_samples.parquet",
        cl_cent=f"{BAGAEV_OUT_DIR}/iatlas_deconv_cluster_centroids_k4.parquet",
        cl_resp=f"{BAGAEV_OUT_DIR}/iatlas_deconv_cluster_response_rates.parquet",
        cl_cont=f"{BAGAEV_OUT_DIR}/iatlas_deconv_clusters_vs_bagaev_contingency.parquet",
        multires=f"{BAGAEV_OUT_DIR}/iatlas_multires_deconv_cluster_benchmark.parquet",
        within=f"{BAGAEV_OUT_DIR}/iatlas_within_subtype_deconv_auc.parquet",
        sc_sum=f"{BAGAEV_OUT_DIR}/sc23_pseudobulk_cohort_summary.parquet",
        sc9_resp=f"{BAGAEV_OUT_DIR}/sc9_milopy_subtype_response_rates.parquet",
        sc9_auc=f"{BAGAEV_OUT_DIR}/sc9_milopy_cohort_bagaev_aucs.parquet",
    output:
        fig1=f"{BAGAEV_RESULTS_DIR}/fig1_iatlas_29fges_portrait_and_clustering.svg",
        fig2=f"{BAGAEV_RESULTS_DIR}/fig2_iatlas_subtype_response_rates.svg",
        fig3=f"{BAGAEV_RESULTS_DIR}/fig3_iatlas_deconvolution_integration.svg",
        fig4=f"{BAGAEV_RESULTS_DIR}/fig4_sc23_pseudobulk_clustering_conservation.svg",
        fig5=f"{BAGAEV_RESULTS_DIR}/fig5_sc9_milopy_response_correlation.svg",
        fig6=f"{BAGAEV_RESULTS_DIR}/fig6_deconv_subtype_correlations.svg",
        fig7=f"{BAGAEV_RESULTS_DIR}/fig7_deconv_sample_clusters_and_response.svg",
    params:
        data_dir=BAGAEV_OUT_DIR,
        results_dir=BAGAEV_RESULTS_DIR,
    threads: 2
    shell:
        "python3 {input.script} --data-dir {params.data_dir} --results-dir {params.results_dir}"
