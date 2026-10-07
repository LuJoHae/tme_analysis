# workflow/rules/reference_sampling_hpo.smk
# Snakemake workflow rules for Single-Cell Reference Sampling and Classifier HPO

HPO_DATA_DIR = config.get("data_dir", DATA_DIR if "DATA_DIR" in locals() else "/storage/halu/data-test")
HPO_OUT_DIR = config.get("hpo_out_dir", f"{HPO_DATA_DIR}/output/reference_sampling_hpo")
HPO_RESULTS_DIR = config.get("hpo_results_dir", f"{HPO_DATA_DIR}/results/reference_sampling_hpo")
HPO_PREPROCESSED_DIR = config.get("hpo_preprocessed_dir", f"{HPO_DATA_DIR}/preprocessed")

HPO_CANCER_TYPES = config.get("hpo_cancer_types", "melanoma")
HPO_N_TRIALS = config.get("hpo_n_trials", 25)
HPO_SEED = config.get("hpo_seed", 42)
HPO_DISCOVERY_COHORTS = config.get(
    "hpo_discovery_cohorts",
    "Hugo,Riaz,Liu",
)
HPO_HELD_OUT_COHORTS = config.get(
    "hpo_held_out_cohorts",
    "Gide,VanAllen,Snyder",
)

rule all_reference_sampling_hpo:
    input:
        f"{HPO_OUT_DIR}/hpo_evaluations.parquet",
        f"{HPO_OUT_DIR}/best_reference/reference_phi.parquet",
        f"{HPO_RESULTS_DIR}/pareto_frontier_loco_vs_collinearity.svg",
        f"{HPO_RESULTS_DIR}/held_out_validation_generalization.svg",

rule run_reference_sampling_hpo:
    input:
        script="scripts/reference_sampling_hpo/run_hpo_pipeline.py",
    output:
        evaluations=f"{HPO_OUT_DIR}/hpo_evaluations.parquet",
        phi=f"{HPO_OUT_DIR}/best_reference/reference_phi.parquet",
        pareto_svg=f"{HPO_RESULTS_DIR}/pareto_frontier_loco_vs_collinearity.svg",
        held_out_svg=f"{HPO_RESULTS_DIR}/held_out_validation_generalization.svg",
    params:
        cancer_types=HPO_CANCER_TYPES,
        n_trials=HPO_N_TRIALS,
        seed=HPO_SEED,
        discovery_cohorts=HPO_DISCOVERY_COHORTS,
        held_out_cohorts=HPO_HELD_OUT_COHORTS,
        preprocessed_dir=HPO_PREPROCESSED_DIR,
        out_dir=HPO_OUT_DIR,
        results_dir=HPO_RESULTS_DIR,
    resources:
        mem_mb=64000,
        cpus=16,
    threads: 16
    shell:
        """
        set -euo pipefail
        ALLOW_LOCAL_HPO=1 python {input.script} \
            --cancer-types {params.cancer_types} \
            --n-trials {params.n_trials} \
            --seed {params.seed} \
            --discovery-cohorts "{params.discovery_cohorts}" \
            --held-out-cohorts "{params.held_out_cohorts}" \
            --preprocessed-dir {params.preprocessed_dir} \
            --out-dir {params.out_dir} \
            --results-dir {params.results_dir} \
            --n-jobs {threads}
        """

rule run_reference_sampling_hpo_benchmark:
    input:
        script="scripts/reference_sampling_hpo/run_hpo_pipeline.py",
    output:
        evaluations=f"{HPO_OUT_DIR}/benchmark/hpo_evaluations.parquet",
        pareto_svg=f"{HPO_RESULTS_DIR}/benchmark/pareto_frontier_loco_vs_collinearity.svg",
    params:
        out_dir=f"{HPO_OUT_DIR}/benchmark",
        results_dir=f"{HPO_RESULTS_DIR}/benchmark",
    threads: 2
    shell:
        """
        set -euo pipefail
        ALLOW_LOCAL_HPO=1 python {input.script} \
            --benchmark-mode \
            --n-trials 3 \
            --out-dir {params.out_dir} \
            --results-dir {params.results_dir} \
            --n-jobs {threads}
        """
