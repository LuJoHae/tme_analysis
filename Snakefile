# Snakefile
# Orchestrates data processing on the remote server

DATA_DIR = "/storage/halu/data"
RESULTS_DIR = f"{DATA_DIR}/results"
GSE_DIR = f"{DATA_DIR}/GSE120575"

# The 'all' rule defines what files should ultimately be generated.
rule all:
    input:
        DATA_DIR + "/GSE120575/gse120575_parsed_cell_ids.parquet",
        DATA_DIR + "/results/whole_dataset/gse120575_umap_plots.svg",
        DATA_DIR + "/results/whole_dataset/gse120575_umap_plots.png",
        DATA_DIR + "/results/whole_dataset/plot_abundance_completed.txt",
        DATA_DIR + "/results/whole_dataset/plot_overlaps_completed.txt",
        DATA_DIR + "/output/whole_dataset/milopy_completed.txt",
        DATA_DIR + "/results/whole_dataset/plot_milopy_completed.txt",
        DATA_DIR + "/results/whole_dataset/plot_nhood_composition_completed.txt",
        DATA_DIR + "/results/tcells/plot_milopy_completed.txt",
        DATA_DIR + "/results/tcells/plot_nhood_composition_completed.txt",
        DATA_DIR + "/results/tcell-CD8/plot_milopy_completed.txt",
        DATA_DIR + "/results/tcell-CD8/plot_nhood_composition_completed.txt",
        DATA_DIR + "/output/tcell-CD8/cd8_responder_dge.csv",
        DATA_DIR + "/results/tcell-CD8/cd8_responder_dge_volcano.svg"

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
        plots_svg = f"{RESULTS_DIR}/whole_dataset/gse120575_umap_plots.svg",
        plots_png = f"{RESULTS_DIR}/whole_dataset/gse120575_umap_plots.png"
    shell:
        """
        mkdir -p {RESULTS_DIR}/whole_dataset
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
        marker=f"{DATA_DIR}/output/whole_dataset/sccoda_completed.txt"
    shell:
        """
        mkdir -p {DATA_DIR}/output/whole_dataset
        python {input.script} \
            --adata {input.adata} \
            --out-dir {DATA_DIR}/output/whole_dataset
        touch {output.marker}
        """

rule plot_abundance:
    input:
        script="scripts/gse120575/plot_cluster_abundance.py",
        marker=f"{DATA_DIR}/output/whole_dataset/sccoda_completed.txt"
    output:
        marker=f"{RESULTS_DIR}/whole_dataset/plot_abundance_completed.txt"
    shell:
        """
        mkdir -p {RESULTS_DIR}/whole_dataset
        python {input.script} \
            --data-dir {DATA_DIR}/output/whole_dataset \
            --out-dir {RESULTS_DIR}/whole_dataset
        touch {output.marker}
        """

rule plot_overlaps:
    input:
        script="scripts/gse120575/plot_cluster_overlaps.py",
        adata=f"{GSE_DIR}/gse120575_processed.h5ad"
    output:
        marker=f"{RESULTS_DIR}/whole_dataset/plot_overlaps_completed.txt"
    shell:
        """
        mkdir -p {RESULTS_DIR}/whole_dataset
        python {input.script} \
            --adata {input.adata} \
            --out-dir {RESULTS_DIR}/whole_dataset
        touch {output.marker}
        """

rule milopy_abundance:
    input:
        script="scripts/gse120575/run_milopy.py",
        adata=f"{GSE_DIR}/gse120575_processed.h5ad"
    output:
        marker=f"{DATA_DIR}/output/whole_dataset/milopy_completed.txt"
    shell:
        """
        mkdir -p {DATA_DIR}/output/whole_dataset
        python {input.script} \
            --adata {input.adata} \
            --out-dir {DATA_DIR}/output/whole_dataset
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

rule plot_milopy:
    input:
        script="scripts/gse120575/plot_milopy.py",
        marker=f"{DATA_DIR}/output/whole_dataset/milopy_completed.txt"
    output:
        marker=f"{RESULTS_DIR}/whole_dataset/plot_milopy_completed.txt"
    shell:
        """
        mkdir -p {RESULTS_DIR}/whole_dataset
        python {input.script} \
            --data-dir {DATA_DIR}/output/whole_dataset \
            --out-dir {RESULTS_DIR}/whole_dataset
        touch {output.marker}
        """

rule plot_nhood_composition:
    input:
        script="scripts/gse120575/plot_nhood_composition.py",
        adata=f"{GSE_DIR}/gse120575_processed.h5ad",
        marker=f"{DATA_DIR}/output/whole_dataset/milopy_completed.txt"
    output:
        marker=f"{RESULTS_DIR}/whole_dataset/plot_nhood_composition_completed.txt"
    shell:
        """
        mkdir -p {RESULTS_DIR}/whole_dataset
        python {input.script} \
            --adata {input.adata} \
            --data-dir {DATA_DIR}/output/whole_dataset \
            --out-dir {RESULTS_DIR}/whole_dataset
        touch {output.marker}
        """

rule check_responder_nhoods:
    input:
        script="scripts/gse120575/check_responder_nhoods.py"
    output:
        marker=f"{RESULTS_DIR}/whole_dataset/check_responder_nhoods_completed.txt",
        plot=f"{RESULTS_DIR}/whole_dataset/neighborhood_fractions_plot.png"
    shell:
        """
        mkdir -p {RESULTS_DIR}/whole_dataset
        python {input.script}
        mv neighborhood_fractions_plot.png {output.plot}
        touch {output.marker}
        """

rule subset_tcells:
    input:
        script="scripts/gse120575/subset_tcells.py",
        adata=f"{GSE_DIR}/gse120575_processed.h5ad"
    output:
        adata=f"{GSE_DIR}/gse120575_tcells_processed.h5ad"
    shell:
        """
        python {input.script} \
            --adata {input.adata} \
            --out-adata {output.adata}
        """

rule milopy_abundance_tcells:
    input:
        script="scripts/gse120575/run_milopy.py",
        adata=f"{GSE_DIR}/gse120575_tcells_processed.h5ad"
    output:
        marker=f"{DATA_DIR}/output/tcells/milopy_completed.txt"
    shell:
        """
        mkdir -p {DATA_DIR}/output/tcells
        python {input.script} \
            --adata {input.adata} \
            --out-dir {DATA_DIR}/output/tcells
        touch {output.marker}
        """

rule plot_milopy_tcells:
    input:
        script="scripts/gse120575/plot_milopy.py",
        marker=f"{DATA_DIR}/output/tcells/milopy_completed.txt"
    output:
        marker=f"{RESULTS_DIR}/tcells/plot_milopy_completed.txt"
    shell:
        """
        mkdir -p {RESULTS_DIR}/tcells
        python {input.script} \
            --data-dir {DATA_DIR}/output/tcells \
            --out-dir {RESULTS_DIR}/tcells
        touch {output.marker}
        """

rule plot_nhood_composition_tcells:
    input:
        script="scripts/gse120575/plot_nhood_composition.py",
        adata=f"{GSE_DIR}/gse120575_tcells_processed.h5ad",
        marker=f"{DATA_DIR}/output/tcells/milopy_completed.txt"
    output:
        marker=f"{RESULTS_DIR}/tcells/plot_nhood_composition_completed.txt"
    shell:
        """
        mkdir -p {RESULTS_DIR}/tcells
        python {input.script} \
            --adata {input.adata} \
            --data-dir {DATA_DIR}/output/tcells \
            --out-dir {RESULTS_DIR}/tcells
        touch {output.marker}
        """

rule subset_tcells_cd8:
    input:
        script="scripts/gse120575/subset_tcell_cd8.py",
        adata=f"{GSE_DIR}/gse120575_tcells_processed.h5ad"
    output:
        adata=f"{GSE_DIR}/gse120575_tcells_cd8_processed.h5ad"
    shell:
        """
        python {input.script} \
            --adata {input.adata} \
            --out-adata {output.adata}
        """

rule milopy_abundance_tcells_cd8:
    input:
        script="scripts/gse120575/run_milopy.py",
        adata=f"{GSE_DIR}/gse120575_tcells_cd8_processed.h5ad"
    output:
        marker=f"{DATA_DIR}/output/tcell-CD8/milopy_completed.txt"
    shell:
        """
        mkdir -p {DATA_DIR}/output/tcell-CD8
        python {input.script} \
            --adata {input.adata} \
            --out-dir {DATA_DIR}/output/tcell-CD8
        touch {output.marker}
        """

rule plot_milopy_tcells_cd8:
    input:
        script="scripts/gse120575/plot_milopy.py",
        marker=f"{DATA_DIR}/output/tcell-CD8/milopy_completed.txt"
    output:
        marker=f"{RESULTS_DIR}/tcell-CD8/plot_milopy_completed.txt"
    shell:
        """
        mkdir -p {RESULTS_DIR}/tcell-CD8
        python {input.script} \
            --data-dir {DATA_DIR}/output/tcell-CD8 \
            --out-dir {RESULTS_DIR}/tcell-CD8
        touch {output.marker}
        """

rule plot_nhood_composition_tcells_cd8:
    input:
        script="scripts/gse120575/plot_nhood_composition.py",
        adata=f"{GSE_DIR}/gse120575_tcells_cd8_processed.h5ad",
        marker=f"{DATA_DIR}/output/tcell-CD8/milopy_completed.txt"
    output:
        marker=f"{RESULTS_DIR}/tcell-CD8/plot_nhood_composition_completed.txt"
    shell:
        """
        mkdir -p {RESULTS_DIR}/tcell-CD8
        python {input.script} \
            --adata {input.adata} \
            --data-dir {DATA_DIR}/output/tcell-CD8 \
            --out-dir {RESULTS_DIR}/tcell-CD8
        touch {output.marker}
        """

rule calc_cd8_dge:
    input:
        script="scripts/gse120575/calc_cd8_dge.py",
        adata=f"{GSE_DIR}/gse120575_tcells_cd8_processed.h5ad"
    output:
        csv=f"{DATA_DIR}/output/tcell-CD8/cd8_responder_dge.csv"
    shell:
        """
        mkdir -p {DATA_DIR}/output/tcell-CD8
        python {input.script} \
            --adata {input.adata} \
            --out-csv {output.csv}
        """

rule plot_cd8_dge:
    input:
        script="scripts/gse120575/plot_cd8_dge.py",
        adata=f"{GSE_DIR}/gse120575_tcells_cd8_processed.h5ad",
        csv=f"{DATA_DIR}/output/tcell-CD8/cd8_responder_dge.csv"
    output:
        plot1=f"{RESULTS_DIR}/tcell-CD8/cd8_responder_dge_volcano.svg",
        plot2=f"{RESULTS_DIR}/tcell-CD8/cd8_responder_dge_dotplot.svg"
    shell:
        """
        mkdir -p {RESULTS_DIR}/tcell-CD8
        python {input.script} \
            --adata {input.adata} \
            --csv {input.csv} \
            --out-dir {RESULTS_DIR}/tcell-CD8
        """

