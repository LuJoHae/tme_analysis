# Snakefile
# Orchestrates data processing on the remote server

DATA_DIR = "/storage/halu/data"
RESULTS_DIR = f"{DATA_DIR}/results"
GSE_DIR = f"{DATA_DIR}/GSE120575"
GSE123139_DIR = f"{DATA_DIR}/GSE123139"
GSE97168_DIR = f"{DATA_DIR}/GSE97168"

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
        DATA_DIR + "/results/tcell-CD8/cd8_responder_dge_volcano.svg",
        DATA_DIR + "/results/tcell-CD4/plot_milopy_completed.txt",
        DATA_DIR + "/results/tcell-CD4/plot_nhood_composition_completed.txt",
        DATA_DIR + "/output/tcell-CD4/cd4_responder_dge.csv",
        DATA_DIR + "/results/tcell-CD4/cd4_responder_dge_volcano.svg",
        DATA_DIR + "/results/myeloid/plot_milopy_completed.txt",
        DATA_DIR + "/output/myeloid/myeloid_responder_dge.csv",
        DATA_DIR + "/results/myeloid/myeloid_responder_dge_volcano.svg",
        DATA_DIR + "/unified_step07_gse120575.h5ad",
        DATA_DIR + "/results/unified_integration/pca_uncorrected_gse_vs_rest.svg",
        DATA_DIR + "/results/unified_integration/umap_uncorrected_gse_vs_rest.svg"

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

rule download_gse123139:
    input:
        script="scripts/gse123139/download_gse123139.py"
    output:
        raw_tar = f"{GSE123139_DIR}/GSE123139_RAW.tar",
        tcr = f"{GSE123139_DIR}/GSE123139_T_cells_tcrb_v2.txt.gz",
        matrix = f"{GSE123139_DIR}/GSE123139_series_matrix.txt.gz"
    params:
        out_dir = GSE123139_DIR
    shell:
        """
        python {input.script} --out-dir {params.out_dir}
        """

rule create_raw_h5ad_gse123139:
    input:
        script="scripts/gse123139/process_gse123139.py",
        raw_tar = f"{GSE123139_DIR}/GSE123139_RAW.tar",
        tcr = f"{GSE123139_DIR}/GSE123139_T_cells_tcrb_v2.txt.gz",
        matrix = f"{GSE123139_DIR}/GSE123139_series_matrix.txt.gz"
    output:
        adata_raw = f"{GSE123139_DIR}/gse123139_raw.h5ad"
    params:
        input_dir = GSE123139_DIR
    shell:
        """
        python {input.script} \
            --input-dir {params.input_dir} \
            --out-h5ad {output.adata_raw}
        """

rule preprocess_gse123139:
    input:
        script="scripts/gse123139/preprocess_gse123139.py",
        adata_raw = f"{GSE123139_DIR}/gse123139_raw.h5ad"
    output:
        adata = f"{GSE123139_DIR}/gse123139_processed.h5ad"
    shell:
        """
        python {input.script} \
            --input-h5ad {input.adata_raw} \
            --out-h5ad {output.adata}
        """

rule plot_gse123139:
    input:
        script="scripts/gse123139/plot_gse123139.py",
        adata = f"{GSE123139_DIR}/gse123139_processed.h5ad"
    output:
        plots_svg = f"{RESULTS_DIR}/gse123139/gse123139_umap_plots.svg",
        plots_png = f"{RESULTS_DIR}/gse123139/gse123139_umap_plots.png"
    shell:
        """
        mkdir -p {RESULTS_DIR}/gse123139
        python {input.script} \
            --adata {input.adata} \
            --out-svg {output.plots_svg} \
            --out-png {output.plots_png}
        """

rule download_gse97168:
    input:
        script="scripts/gse97168/download_gse97168.py"
    output:
        umitab = f"{GSE97168_DIR}/GSE97168_umitab.txt.gz",
        meta = f"{GSE97168_DIR}/GSE97168_metadata.txt.gz"
    params:
        out_dir = GSE97168_DIR
    shell:
        """
        python {input.script} --out-dir {params.out_dir}
        """

rule create_raw_h5ad_gse97168:
    input:
        script="scripts/gse97168/process_gse97168.py",
        umitab = f"{GSE97168_DIR}/GSE97168_umitab.txt.gz",
        meta = f"{GSE97168_DIR}/GSE97168_metadata.txt.gz"
    output:
        adata_raw = f"{GSE97168_DIR}/gse97168_raw.h5ad"
    params:
        input_dir = GSE97168_DIR
    shell:
        """
        python {input.script} \
            --input-dir {params.input_dir} \
            --out-h5ad {output.adata_raw}
        """

rule preprocess_gse97168:
    input:
        script="scripts/gse97168/preprocess_gse97168.py",
        adata_raw = f"{GSE97168_DIR}/gse97168_raw.h5ad"
    output:
        adata = f"{GSE97168_DIR}/gse97168_processed.h5ad"
    shell:
        """
        python {input.script} \
            --input-h5ad {input.adata_raw} \
            --out-h5ad {output.adata}
        """

rule plot_gse97168:
    input:
        script="scripts/gse97168/plot_gse97168.py",
        adata = f"{GSE97168_DIR}/gse97168_processed.h5ad"
    output:
        plots_svg = f"{RESULTS_DIR}/gse97168/gse97168_umap_plots.svg",
        plots_png = f"{RESULTS_DIR}/gse97168/gse97168_umap_plots.png"
    shell:
        """
        mkdir -p {RESULTS_DIR}/gse97168
        python {input.script} \
            --adata {input.adata} \
            --out-svg {output.plots_svg} \
            --out-png {output.plots_png}
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
        script="scripts/gse120575/subset_lineage.py",
        adata=f"{GSE_DIR}/gse120575_processed.h5ad"
    output:
        adata=f"{GSE_DIR}/gse120575_tcells_processed.h5ad"
    shell:
        """
        python {input.script} \
            --adata {input.adata} \
            --out-adata {output.adata} \
            --marker CD3D
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
        script="scripts/gse120575/subset_tcell_by_marker.py",
        adata=f"{GSE_DIR}/gse120575_tcells_processed.h5ad"
    output:
        adata=f"{GSE_DIR}/gse120575_tcells_cd8_processed.h5ad"
    shell:
        """
        python {input.script} \
            --adata {input.adata} \
            --out-adata {output.adata} \
            --marker CD8A
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
        script="scripts/gse120575/calc_responder_dge.py",
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
        script="scripts/gse120575/plot_responder_dge.py",
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
            --out-dir {RESULTS_DIR}/tcell-CD8 \
            --prefix cd8 \
            --title 'CD8+ T cells'
        """

rule subset_tcells_cd4:
    input:
        script="scripts/gse120575/subset_tcell_by_marker.py",
        adata=f"{GSE_DIR}/gse120575_tcells_processed.h5ad"
    output:
        adata=f"{GSE_DIR}/gse120575_tcells_cd4_processed.h5ad"
    shell:
        """
        python {input.script} \
            --adata {input.adata} \
            --out-adata {output.adata} \
            --marker CD4
        """

rule milopy_abundance_tcells_cd4:
    input:
        script="scripts/gse120575/run_milopy.py",
        adata=f"{GSE_DIR}/gse120575_tcells_cd4_processed.h5ad"
    output:
        marker=f"{DATA_DIR}/output/tcell-CD4/milopy_completed.txt"
    shell:
        """
        mkdir -p {DATA_DIR}/output/tcell-CD4
        python {input.script} \
            --adata {input.adata} \
            --out-dir {DATA_DIR}/output/tcell-CD4
        touch {output.marker}
        """

rule plot_milopy_tcells_cd4:
    input:
        script="scripts/gse120575/plot_milopy.py",
        marker=f"{DATA_DIR}/output/tcell-CD4/milopy_completed.txt"
    output:
        marker=f"{RESULTS_DIR}/tcell-CD4/plot_milopy_completed.txt"
    shell:
        """
        mkdir -p {RESULTS_DIR}/tcell-CD4
        python {input.script} \
            --data-dir {DATA_DIR}/output/tcell-CD4 \
            --out-dir {RESULTS_DIR}/tcell-CD4
        touch {output.marker}
        """

rule plot_nhood_composition_tcells_cd4:
    input:
        script="scripts/gse120575/plot_nhood_composition.py",
        adata=f"{GSE_DIR}/gse120575_tcells_cd4_processed.h5ad",
        marker=f"{DATA_DIR}/output/tcell-CD4/milopy_completed.txt"
    output:
        marker=f"{RESULTS_DIR}/tcell-CD4/plot_nhood_composition_completed.txt"
    shell:
        """
        mkdir -p {RESULTS_DIR}/tcell-CD4
        python {input.script} \
            --adata {input.adata} \
            --data-dir {DATA_DIR}/output/tcell-CD4 \
            --out-dir {RESULTS_DIR}/tcell-CD4
        touch {output.marker}
        """

rule calc_cd4_dge:
    input:
        script="scripts/gse120575/calc_responder_dge.py",
        adata=f"{GSE_DIR}/gse120575_tcells_cd4_processed.h5ad"
    output:
        csv=f"{DATA_DIR}/output/tcell-CD4/cd4_responder_dge.csv"
    shell:
        """
        mkdir -p {DATA_DIR}/output/tcell-CD4
        python {input.script} \
            --adata {input.adata} \
            --out-csv {output.csv}
        """

rule plot_cd4_dge:
    input:
        script="scripts/gse120575/plot_responder_dge.py",
        adata=f"{GSE_DIR}/gse120575_tcells_cd4_processed.h5ad",
        csv=f"{DATA_DIR}/output/tcell-CD4/cd4_responder_dge.csv"
    output:
        plot1=f"{RESULTS_DIR}/tcell-CD4/cd4_responder_dge_volcano.svg",
        plot2=f"{RESULTS_DIR}/tcell-CD4/cd4_responder_dge_dotplot.svg"
    shell:
        """
        mkdir -p {RESULTS_DIR}/tcell-CD4
        python {input.script} \
            --adata {input.adata} \
            --csv {input.csv} \
            --out-dir {RESULTS_DIR}/tcell-CD4 \
            --prefix cd4 \
            --title 'CD4+ T cells'
        """

rule subset_myeloid:
    input:
        script="scripts/gse120575/subset_lineage.py",
        adata=f"{GSE_DIR}/gse120575_processed.h5ad"
    output:
        adata=f"{GSE_DIR}/gse120575_myeloid_processed.h5ad"
    shell:
        """
        python {input.script} \
            --adata {input.adata} \
            --out-adata {output.adata} \
            --marker ITGAM \
            --exclude CD3D
        """

rule milopy_abundance_myeloid:
    input:
        script="scripts/gse120575/run_milopy.py",
        adata=f"{GSE_DIR}/gse120575_myeloid_processed.h5ad"
    output:
        marker=f"{DATA_DIR}/output/myeloid/milopy_completed.txt"
    shell:
        """
        mkdir -p {DATA_DIR}/output/myeloid
        python {input.script} \
            --adata {input.adata} \
            --out-dir {DATA_DIR}/output/myeloid \
            --fdr 0.2
        touch {output.marker}
        """

rule plot_milopy_myeloid:
    input:
        script="scripts/gse120575/plot_milopy.py",
        marker=f"{DATA_DIR}/output/myeloid/milopy_completed.txt"
    output:
        marker=f"{RESULTS_DIR}/myeloid/plot_milopy_completed.txt"
    shell:
        """
        mkdir -p {RESULTS_DIR}/myeloid
        python {input.script} \
            --data-dir {DATA_DIR}/output/myeloid \
            --out-dir {RESULTS_DIR}/myeloid \
            --fdr 0.2
        touch {output.marker}
        """

rule plot_nhood_composition_myeloid:
    input:
        script="scripts/gse120575/plot_nhood_composition.py",
        adata=f"{GSE_DIR}/gse120575_myeloid_processed.h5ad",
        marker=f"{DATA_DIR}/output/myeloid/milopy_completed.txt"
    output:
        marker=f"{RESULTS_DIR}/myeloid/plot_nhood_composition_completed.txt"
    shell:
        """
        mkdir -p {RESULTS_DIR}/myeloid
        python {input.script} \
            --adata {input.adata} \
            --data-dir {DATA_DIR}/output/myeloid \
            --out-dir {RESULTS_DIR}/myeloid \
            --fdr 0.2
        touch {output.marker}
        """

rule calc_myeloid_dge:
    input:
        script="scripts/gse120575/calc_responder_dge.py",
        adata=f"{GSE_DIR}/gse120575_myeloid_processed.h5ad"
    output:
        csv=f"{DATA_DIR}/output/myeloid/myeloid_responder_dge.csv"
    shell:
        """
        mkdir -p {DATA_DIR}/output/myeloid
        python {input.script} \
            --adata {input.adata} \
            --out-csv {output.csv}
        """

rule plot_myeloid_dge:
    input:
        script="scripts/gse120575/plot_responder_dge.py",
        adata=f"{GSE_DIR}/gse120575_myeloid_processed.h5ad",
        csv=f"{DATA_DIR}/output/myeloid/myeloid_responder_dge.csv"
    output:
        plot1=f"{RESULTS_DIR}/myeloid/myeloid_responder_dge_volcano.svg",
        plot2=f"{RESULTS_DIR}/myeloid/myeloid_responder_dge_dotplot.svg"
    shell:
        """
        mkdir -p {RESULTS_DIR}/myeloid
        python {input.script} \
            --adata {input.adata} \
            --csv {input.csv} \
            --out-dir {RESULTS_DIR}/myeloid \
            --prefix myeloid \
            --title 'Myeloid Cells'
        """

rule unify_step07_gse120575:
    input:
        script="scripts/unify_step07_gse120575.py",
        tpm=f"{GSE_DIR}/gse120575_tpm.parquet",
        tpm_meta=f"{GSE_DIR}/gse120575_tpm_cell_metadata.parquet"
    output:
        h5ad=f"{DATA_DIR}/unified_step07_gse120575.h5ad"
    shell:
        """
        python {input.script} \
            --lair-dir /storage/halu/lair \
            --gse120575-dir {GSE_DIR} \
            --output-path {output.h5ad}
        """

rule plot_unified_embeddings:
    input:
        script="scripts/plot_unified_embeddings.py",
        h5ad=f"{DATA_DIR}/unified_step07_gse120575.h5ad"
    output:
        pca_gse=f"{RESULTS_DIR}/unified_integration/pca_uncorrected_gse_vs_rest.svg",
        umap_gse=f"{RESULTS_DIR}/unified_integration/umap_uncorrected_gse_vs_rest.svg"
    shell:
        """
        mkdir -p {RESULTS_DIR}/unified_integration
        python {input.script} \
            --input-h5ad {input.h5ad} \
            --output-dir {RESULTS_DIR}/unified_integration
        """

rule download_geo_all:
    input:
        script="scripts/01_download_geo_data.py"
    output:
        marker=f"{DATA_DIR}/raw_geo/download_completed.txt"
    shell:
        """
        python {input.script} --out-dir {DATA_DIR}/raw_geo
        touch {output.marker}
        """

rule process_raw_to_sparse_h5:
    input:
        script="scripts/02_process_raw_to_sparse_h5.py",
        marker=f"{DATA_DIR}/raw_geo/download_completed.txt"
    output:
        marker=f"{DATA_DIR}/sparse_h5/conversion_completed.txt"
    shell:
        """
        python {input.script} --raw-dir {DATA_DIR}/raw_geo --out-dir {DATA_DIR}/sparse_h5
        touch {output.marker}
        """

rule preprocess_datasets_to_h5:
    input:
        script="scripts/03_preprocess_datasets_to_h5.py",
        marker=f"{DATA_DIR}/sparse_h5/conversion_completed.txt"
    output:
        marker=f"{DATA_DIR}/processed_h5/preprocessing_completed.txt"
    shell:
        """
        python {input.script} --input-dir {DATA_DIR}/sparse_h5 --out-dir {DATA_DIR}/processed_h5
        touch {output.marker}
        """

