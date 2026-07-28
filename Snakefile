# Snakefile
# Orchestrates data processing on the remote server

DATA_DIR = "/storage/halu/data"
RESULTS_DIR = f"{DATA_DIR}/results"
GSE_DIR = f"{DATA_DIR}/GSE120575"

# The 'all' rule defines what files should ultimately be generated.
rule all:
    input:
        f"{RESULTS_DIR}/gse120575_summary.csv",
        f"{RESULTS_DIR}/gse120575_plots.html",
        f"{RESULTS_DIR}/gse120575_plots.svg",
        f"{GSE_DIR}/gse120575_parsed_cell_ids.parquet"

rule download_gse120575:
    input:
        script="scripts/download_gse120575.py"
    output:
        tpm = f"{GSE_DIR}/gse120575_tpm.parquet",
        tpm_meta = f"{GSE_DIR}/gse120575_tpm_cell_metadata.parquet",
        meta = f"{GSE_DIR}/gse120575_meta.parquet"
    params:
        out_dir = GSE_DIR
    shell:
        """
        python scripts/download_gse120575.py --out-dir {params.out_dir}
        """

rule analyze_gse120575:
    input:
        script="scripts/analyze_gse120575.py",
        tpm = f"{GSE_DIR}/gse120575_tpm.parquet",
        meta = f"{GSE_DIR}/gse120575_meta.parquet"
    output:
        summary = f"{RESULTS_DIR}/gse120575_summary.csv",
        plots = f"{RESULTS_DIR}/gse120575_plots.html",
        plots_svg = f"{RESULTS_DIR}/gse120575_plots.svg"
    shell:
        """
        python scripts/analyze_gse120575.py \
            --tpm {input.tpm} \
            --meta {input.meta} \
            --out-csv {output.summary} \
            --out-html {output.plots} \
            --out-svg {output.plots_svg}
        """

rule extract_metadata_gse120575:
    input:
        script="scripts/extract_cell_metadata.py",
        meta=f"{GSE_DIR}/gse120575_meta.parquet"
    output:
        parsed_meta=f"{GSE_DIR}/gse120575_parsed_cell_ids.parquet"
    shell:
        """
        python {input.script} --meta {input.meta} --out {output.parsed_meta}
        """
