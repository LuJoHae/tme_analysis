# Snakefile
# Orchestrates data processing on the remote server

DATA_DIR = "/storage/halu/data"
RESULTS_DIR = f"{DATA_DIR}/results"
GSE_DIR = f"{DATA_DIR}/GSE120575"

# The 'all' rule defines what files should ultimately be generated.
rule all:
    input:
        DATA_DIR + "/GSE120575/gse120575_parsed_cell_ids.parquet",
        DATA_DIR + "/results/gse120575_umap_plots.svg",
        DATA_DIR + "/results/gse120575_umap_plots.png",
        DATA_DIR + "/results/plot_abundance_completed.txt"

rule download_gse120575:
    input:
        script="scripts/gse120575/download_gse120575.py"
    output:
        tpm = f"{GSE_DIR}/gse120575_tpm.parquet",
        tpm_meta = f"{GSE_DIR}/gse120575_tpm_cell_metadata.parquet",
        meta = f"{GSE_DIR}/gse120575_meta.parquet"
    params:
        out_dir = GSE_DIR
    shell:
        """
        python scripts/gse120575/download_gse120575.py --out-dir {params.out_dir}
        """

rule preprocess_gse120575:
    input:
        script="scripts/gse120575/preprocess_gse120575.py",
        tpm = f"{GSE_DIR}/gse120575_tpm.parquet",
        meta = f"{GSE_DIR}/gse120575_parsed_cell_ids.parquet"
    output:
        adata = f"{GSE_DIR}/gse120575_processed.h5ad"
    shell:
        """
        python {input.script} \
            --tpm {input.tpm} \
            --meta {input.meta} \
            --out-h5ad {output.adata}
        """

rule plot_gse120575:
    input:
        script="scripts/gse120575/plot_gse120575.py",
        adata = f"{GSE_DIR}/gse120575_processed.h5ad"
    output:
        plots_svg = f"{RESULTS_DIR}/gse120575_umap_plots.svg",
        plots_png = f"{RESULTS_DIR}/gse120575_umap_plots.png"
    shell:
        """
        python {input.script} \
            --adata {input.adata} \
            --out-svg {output.plots_svg} \
            --out-png {output.plots_png}
        """

rule cluster_abundance:
    input:
        script="scripts/gse120575/cluster_abundance.py",
        adata=f"{GSE_DIR}/gse120575_processed.h5ad"
    output:
        marker=f"{RESULTS_DIR}/sccoda_completed.txt"
    shell:
        """
        python {input.script} \
            --adata {input.adata} \
            --out-dir {RESULTS_DIR}
        touch {output.marker}
        """

rule plot_abundance:
    input:
        script="scripts/gse120575/plot_cluster_abundance.py",
        marker=f"{RESULTS_DIR}/sccoda_completed.txt"
    output:
        marker=f"{RESULTS_DIR}/plot_abundance_completed.txt"
    shell:
        """
        python {input.script} \
            --data-dir {RESULTS_DIR} \
            --out-dir {RESULTS_DIR}
        touch {output.marker}
        """

rule extract_metadata_gse120575:
    input:
        script="scripts/gse120575/extract_cell_metadata.py",
        meta=f"{GSE_DIR}/gse120575_meta.parquet"
    output:
        parsed_meta=f"{GSE_DIR}/gse120575_parsed_cell_ids.parquet"
    shell:
        """
        python {input.script} --meta {input.meta} --out {output.parsed_meta}
        """
