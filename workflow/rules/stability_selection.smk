# workflow/rules/stability_selection.smk
# Snakemake workflow rules for Multi-Cohort Stability Selection Benchmark & Publication Figures

STABILITY_OUT_DIR = config.get("stability_out_dir", "output/stability_selection_fitters_benchmark")
STABILITY_COHORTS = config.get(
    "stability_cohorts",
    "pancancer melanoma rcc Rosenberg-iAtlas McDermott-iAtlas Riaz-iAtlas Gide-iAtlas Liu-iAtlas Hugo-iAtlas"
)
STABILITY_CUTOFF = config.get("stability_cutoff", 0.75)
STABILITY_Q_BUDGET = config.get("stability_q_budget", 20.0)
STABILITY_B = config.get("stability_B", 50)
STABILITY_SEED = config.get("stability_seed", 42)

rule all_stability_selection:
    input:
        f"{STABILITY_OUT_DIR}/stability_scores.parquet",
        f"{STABILITY_OUT_DIR}/stability_paths.parquet",
        f"{STABILITY_OUT_DIR}/cohort_fitter_summary.parquet",
        f"{STABILITY_OUT_DIR}/figures/fitter_gene_overlap_nature.svg",
        f"{STABILITY_OUT_DIR}/figures/fitter_gene_overlap_nature.png"

rule run_multi_cohort_stability_benchmark:
    input:
        script="scripts/stability_selection_iatlas/run_stability_selection_iatlas.py"
    output:
        scores=f"{STABILITY_OUT_DIR}/stability_scores.parquet",
        paths=f"{STABILITY_OUT_DIR}/stability_paths.parquet",
        summary=f"{STABILITY_OUT_DIR}/cohort_fitter_summary.parquet"
    params:
        cohorts=STABILITY_COHORTS,
        cutoff=STABILITY_CUTOFF,
        q_budget=STABILITY_Q_BUDGET,
        b=STABILITY_B,
        seed=STABILITY_SEED,
        out_dir=STABILITY_OUT_DIR
    shell:
        """
        python {input.script} \
            --cohorts {params.cohorts} \
            --cutoff {params.cutoff} \
            --q-budget {params.q_budget} \
            --B {params.b} \
            --seed {params.seed} \
            --output-dir {params.out_dir}
        """

rule plot_stability_selection_fitter_overlap:
    input:
        script="scripts/stability_selection_iatlas/plotting/plot_fitter_gene_overlap.py",
        scores=f"{STABILITY_OUT_DIR}/stability_scores.parquet"
    output:
        svg=f"{STABILITY_OUT_DIR}/figures/fitter_gene_overlap_nature.svg",
        png=f"{STABILITY_OUT_DIR}/figures/fitter_gene_overlap_nature.png"
    shell:
        """
        python {input.script} \
            --scores-path {input.scores} \
            --output-svg {output.svg} \
            --output-png {output.png}
        """
