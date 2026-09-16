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

DATASETS = [
    "AziziSingleCellMapDiverse2018Adata",
    "BeckerSinglecellAnalysesDefine2022Adata",
    "BiermannDissectingTreatmentnaiveEcosystem2022Adata",
    "BorcherdingMappingImmuneEnvironment2021Adata",
    "ChengPancancerSinglecellTranscriptional2021Adata",
    "DuranteSinglecellAnalysisReveals2020Adata",
    "JerbyArnonCancerCellProgram2018Adata",
    "KhaliqRefiningColorectalCancer2022Adata",
    "KimSinglecellRNASequencing2020Adata",
    "LeaderSinglecellAnalysisHuman2021Adata",
    "LuSinglecellAtlasMulticellular2022Adata",
    "PelkaSpatiallyOrganizedMulticellular2021Adata",
    "PuSinglecellTranscriptomicAnalysis2021Adata",
    "QianPancancerBlueprintHeterogeneous2020aAdata",
    "SharmaOncofetalReprogrammingEndothelial2020Adata",
    "SadeFeldmanDefiningTCell2018Adata",
    "YostClonalReplacementTumor2019Adata",
    "ZhengLandscapeInfiltratingTCells2017Adata",
]

rule download_dataset:
    input:
        script="scripts/01_download_datasets.py"
    output:
        marker=f"{DATA_DIR}/raw_geo/{{dataset}}/download_completed.txt"
    shell:
        """
        python {input.script} --dataset {wildcards.dataset} --out-dir {DATA_DIR}/raw_geo
        touch {output.marker}
        """

rule convert_dataset_to_sparse_h5:
    input:
        script="scripts/02_convert_to_sparse_h5.py",
        marker=f"{DATA_DIR}/raw_geo/{{dataset}}/download_completed.txt"
    output:
        h5ad=f"{DATA_DIR}/sparse_h5/{{dataset}}_sparse.h5ad"
    shell:
        """
        python {input.script} --dataset {wildcards.dataset} --raw-dir {DATA_DIR}/raw_geo --out-h5 {output.h5ad}
        """

rule preprocess_dataset_to_h5:
    input:
        script="scripts/03_preprocess_datasets.py",
        h5ad=f"{DATA_DIR}/sparse_h5/{{dataset}}_sparse.h5ad"
    output:
        h5ad=f"{DATA_DIR}/processed_h5/{{dataset}}_processed.h5ad"
    shell:
        """
        python {input.script} --input-h5 {input.h5ad} --out-h5 {output.h5ad}
        """

rule plot_dataset_embedding:
    input:
        script="scripts/04_plot_dataset_embeddings.py",
        h5ad=f"{DATA_DIR}/processed_h5/{{dataset}}_processed.h5ad"
    output:
        marker=f"{RESULTS_DIR}/results_embeddings/{{dataset}}_completed.txt"
    shell:
        """
        mkdir -p {RESULTS_DIR}/results_embeddings
        python {input.script} --input-h5 {input.h5ad} --out-dir {RESULTS_DIR}/results_embeddings
        touch {output.marker}
        """

rule plot_all_dataset_embeddings:
    input:
        expand(RESULTS_DIR + "/results_embeddings/{dataset}_completed.txt", dataset=DATASETS)


# ==============================================================================
# Sade-Feldman Deconvolution & Milopy Validation Pipeline
# ==============================================================================
SADE_VALIDATION_OUT = f"{DATA_DIR}/output/sade_feldman_deconv_validation"
SADE_VALIDATION_RES = f"{RESULTS_DIR}/sade_feldman_deconv_validation"

rule build_sade_feldman_reference:
    input:
        script="scripts/sade_feldman_deconv_validation/01_build_reference.py",
        adata=f"{GSE_DIR}/gse120575_processed.h5ad",
        tpm=f"{GSE_DIR}/gse120575_tpm.parquet",
    output:
        phi=f"{SADE_VALIDATION_OUT}/reference_phi.parquet",
        markers=f"{SADE_VALIDATION_OUT}/reference_marker_genes.parquet",
    params:
        out_dir=SADE_VALIDATION_OUT,
    shell:
        """
        python {input.script} --adata {input.adata} --tpm {input.tpm} --out-dir {params.out_dir}
        """

rule deconvolute_iatlas_sade_feldman:
    input:
        script="scripts/sade_feldman_deconv_validation/02_deconvolute_iatlas.py",
        phi=f"{SADE_VALIDATION_OUT}/reference_phi.parquet",
    output:
        fracs=f"{SADE_VALIDATION_OUT}/deconv_fractions.parquet",
    params:
        out_dir=SADE_VALIDATION_OUT,
        lair_dir="/storage/halu/lair",
    shell:
        """
        python {input.script} --reference {input.phi} --lair-dir {params.lair_dir} --out-dir {params.out_dir}
        """

rule logistic_regression_sade_feldman:
    input:
        script="scripts/sade_feldman_deconv_validation/03_logistic_regression.py",
        fracs=f"{SADE_VALIDATION_OUT}/deconv_fractions.parquet",
    output:
        results=f"{SADE_VALIDATION_OUT}/logistic_regression_results.parquet",
    params:
        out_dir=SADE_VALIDATION_OUT,
        lair_dir="/storage/halu/lair",
    shell:
        """
        python {input.script} --fractions {input.fracs} --lair-dir {params.lair_dir} --out-dir {params.out_dir}
        """

rule analyze_milopy_sade_feldman:
    input:
        script="scripts/sade_feldman_deconv_validation/04_analyze_milopy.py",
        adata=f"{GSE_DIR}/gse120575_processed.h5ad",
    output:
        states=f"{SADE_VALIDATION_OUT}/milopy_cell_state_da.parquet",
        all_states=f"{SADE_VALIDATION_OUT}/milopy_cell_state_da_all.parquet",
        pre_states=f"{SADE_VALIDATION_OUT}/milopy_cell_state_da_Pre.parquet",
        post_states=f"{SADE_VALIDATION_OUT}/milopy_cell_state_da_Post.parquet",
        cells=f"{SADE_VALIDATION_OUT}/milopy_cell_level_scores.parquet",
    params:
        out_dir=SADE_VALIDATION_OUT,
        milo_dir=f"{DATA_DIR}/output/whole_dataset",
    shell:
        """
        python {input.script} --adata {input.adata} --milo-dir {params.milo_dir} --out-dir {params.out_dir}
        """

rule build_integrated_reference:
    input:
        script="scripts/sade_feldman_deconv_validation/01b_build_integrated_reference.py",
        adata=f"{GSE_DIR}/gse120575_processed.h5ad",
    output:
        phi=f"{SADE_VALIDATION_OUT}/integrated_reference_phi.parquet",
        markers=f"{SADE_VALIDATION_OUT}/integrated_marker_genes.parquet",
        meta=f"{SADE_VALIDATION_OUT}/integrated_cell_metadata.parquet",
        metrics=f"{SADE_VALIDATION_OUT}/integration_quality_metrics.parquet",
    params:
        out_dir=SADE_VALIDATION_OUT,
        lair_dir="/storage/halu/lair",
    shell:
        """
        python {input.script} --adata-sf {input.adata} --lair-dir {params.lair_dir} --out-dir {params.out_dir}
        """

rule concordance_sade_feldman:
    input:
        script="scripts/sade_feldman_deconv_validation/05_compare_concordance.py",
        logistic=f"{SADE_VALIDATION_OUT}/logistic_regression_results.parquet",
        milo_all=f"{SADE_VALIDATION_OUT}/milopy_cell_state_da_all.parquet",
    output:
        metrics=f"{SADE_VALIDATION_OUT}/concordance_metrics.parquet",
        summary=f"{SADE_VALIDATION_OUT}/concordance_summary.parquet",
        cohorts=f"{SADE_VALIDATION_OUT}/cohort_level_concordance_summary.parquet",
        timepoints=f"{SADE_VALIDATION_OUT}/timepoint_concordance_summary.parquet",
        full_metrics=f"{SADE_VALIDATION_OUT}/concordance_metrics_full.parquet",
    params:
        out_dir=SADE_VALIDATION_OUT,
    shell:
        """
        python {input.script} --logistic-results {input.logistic} --milo-dir {params.out_dir} --out-dir {params.out_dir}
        """

rule plot_sade_feldman_validation:
    input:
        script="scripts/sade_feldman_deconv_validation/06_plot_figures.py",
        markers=f"{SADE_VALIDATION_OUT}/reference_marker_genes.parquet",
        fracs=f"{SADE_VALIDATION_OUT}/deconv_fractions.parquet",
        logistic=f"{SADE_VALIDATION_OUT}/logistic_regression_results.parquet",
        milo_all=f"{SADE_VALIDATION_OUT}/milopy_cell_state_da_all.parquet",
        cohorts=f"{SADE_VALIDATION_OUT}/cohort_level_concordance_summary.parquet",
        metrics_full=f"{SADE_VALIDATION_OUT}/concordance_metrics_full.parquet",
        cells=f"{SADE_VALIDATION_OUT}/milopy_cell_level_scores.parquet",
        integ_meta=f"{SADE_VALIDATION_OUT}/integrated_cell_metadata.parquet",
    output:
        fig1=f"{SADE_VALIDATION_RES}/step01_reference_marker_heatmap.svg",
        fig1b=f"{SADE_VALIDATION_RES}/step01b_integrated_umap_batch_correction.svg",
        fig2=f"{SADE_VALIDATION_RES}/step02_deconv_fractions_distribution.svg",
        fig2b=f"{SADE_VALIDATION_RES}/step02b_cohort_deconv_fractions_grid.svg",
        fig2c=f"{SADE_VALIDATION_RES}/step02c_cohort_cell_state_stacked_bars.svg",
        fig3a=f"{SADE_VALIDATION_RES}/step03_logistic_regression_volcano.svg",
        fig3b=f"{SADE_VALIDATION_RES}/step03_logistic_regression_forest.svg",
        fig4a=f"{SADE_VALIDATION_RES}/step04_milopy_nhood_volcano.svg",
        fig4b=f"{SADE_VALIDATION_RES}/step04_milopy_cell_state_da.svg",
        fig4c=f"{SADE_VALIDATION_RES}/step04b_milopy_pre_vs_post.svg",
        fig5=f"{SADE_VALIDATION_RES}/step05_concordance_scatter.svg",
        fig5b=f"{SADE_VALIDATION_RES}/step05b_cohort_concordance_comparison.svg",
        fig5c=f"{SADE_VALIDATION_RES}/step05b_individual_cohort_scatters.svg",
        fig6=f"{SADE_VALIDATION_RES}/step06_dual_umap_validation.svg",
        fig6b=f"{SADE_VALIDATION_RES}/step06b_stratified_umap_validation.svg",
    params:
        data_dir=SADE_VALIDATION_OUT,
        results_dir=SADE_VALIDATION_RES,
    shell:
        """
        python {input.script} --data-dir {params.data_dir} --results-dir {params.results_dir}
        """

rule run_extended_sade_feldman_pipeline:
    input:
        f"{SADE_VALIDATION_RES}/step01_reference_marker_heatmap.svg",
        f"{SADE_VALIDATION_RES}/step01b_integrated_umap_batch_correction.svg",
        f"{SADE_VALIDATION_RES}/step02_deconv_fractions_distribution.svg",
        f"{SADE_VALIDATION_RES}/step02b_cohort_deconv_fractions_grid.svg",
        f"{SADE_VALIDATION_RES}/step02c_cohort_cell_state_stacked_bars.svg",
        f"{SADE_VALIDATION_RES}/step03_logistic_regression_volcano.svg",
        f"{SADE_VALIDATION_RES}/step03_logistic_regression_forest.svg",
        f"{SADE_VALIDATION_RES}/step04_milopy_nhood_volcano.svg",
        f"{SADE_VALIDATION_RES}/step04_milopy_cell_state_da.svg",
        f"{SADE_VALIDATION_RES}/step04b_milopy_pre_vs_post.svg",
        f"{SADE_VALIDATION_RES}/step05_concordance_scatter.svg",
        f"{SADE_VALIDATION_RES}/step05b_cohort_concordance_comparison.svg",
        f"{SADE_VALIDATION_RES}/step05b_individual_cohort_scatters.svg",
        f"{SADE_VALIDATION_RES}/step06_dual_umap_validation.svg",
        f"{SADE_VALIDATION_RES}/step06b_stratified_umap_validation.svg",

rule run_sade_feldman_pipeline:
    input:
        rules.run_extended_sade_feldman_pipeline.input,




