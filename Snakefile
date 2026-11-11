# Snakefile
# Orchestrates data processing on the remote CSCS Alps Eiger cluster via Slurm

import os

DATA_DIR = os.environ.get("TME_DATA_DIR", config.get("data_dir", "/ritom/scratch/cscs/lhaeuser/data-test"))
if str(DATA_DIR).startswith("/capstor/store"):
    raise RuntimeError(
        f"Direct compute read/write to persistent $STORE ('{DATA_DIR}') is strictly prohibited. "
        "Stage data to $SCRATCH (/ritom/scratch/cscs/lhaeuser/data-test) first via 'make stage-in'."
    )

LAIR_DIR = os.environ.get("TME_LAIR_DIR", f"{DATA_DIR}/lair")
RESULTS_DIR = f"{DATA_DIR}/results"
GSE_DIR = f"{DATA_DIR}/GSE120575"
GSE123139_DIR = f"{DATA_DIR}/GSE123139"
GSE97168_DIR = f"{DATA_DIR}/GSE97168"
MANUSCRIPT_OUT_DIR = os.environ.get(
    "TME_PIPELINE_OUTPUT_ROOT", config.get("manuscript_out_dir", "output/manuscript_pipeline")
)

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
        if [ "${{TME_DISABLE_DOWNLOADS:-0}}" = "1" ] || [ -d /ritom/scratch ]; then
            echo "ERROR: Internet dataset downloads are strictly forbidden on CSCS Eiger. Download on olm and transfer via Globus." >&2
            exit 1
        fi
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
        if [ "${{TME_DISABLE_DOWNLOADS:-0}}" = "1" ] || [ -d /ritom/scratch ]; then
            echo "ERROR: Internet dataset downloads are strictly forbidden on CSCS Eiger. Download on olm and transfer via Globus." >&2
            exit 1
        fi
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
        if [ "${{TME_DISABLE_DOWNLOADS:-0}}" = "1" ] || [ -d /ritom/scratch ]; then
            echo "ERROR: Internet dataset downloads are strictly forbidden on CSCS Eiger. Download on olm and transfer via Globus." >&2
            exit 1
        fi
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
        script="scripts/ingestion/unify_step07_gse120575.py",
        tpm=f"{GSE_DIR}/gse120575_tpm.parquet",
        tpm_meta=f"{GSE_DIR}/gse120575_tpm_cell_metadata.parquet"
    output:
        h5ad=f"{DATA_DIR}/unified_step07_gse120575.h5ad"
    params:
        lair_dir=LAIR_DIR
    shell:
        """
        python {input.script} \
            --lair-dir {params.lair_dir} \
            --gse120575-dir {GSE_DIR} \
            --output-path {output.h5ad}
        """

rule plot_unified_embeddings:
    input:
        script="scripts/ingestion/plot_unified_embeddings.py",
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
        script="scripts/ingestion/01_download_datasets.py"
    output:
        marker=f"{DATA_DIR}/raw_geo/{{dataset}}/download_completed.txt"
    shell:
        """
        if [ "${{TME_DISABLE_DOWNLOADS:-0}}" = "1" ] || [ -d /ritom/scratch ]; then
            echo "ERROR: Internet dataset downloads are strictly forbidden on CSCS Eiger. Download on olm and transfer via Globus." >&2
            exit 1
        fi
        python {input.script} --dataset {wildcards.dataset} --out-dir {DATA_DIR}/raw_geo
        touch {output.marker}
        """

rule convert_dataset_to_sparse_h5:
    input:
        script="scripts/ingestion/02_convert_to_sparse_h5.py",
        marker=f"{DATA_DIR}/raw_geo/{{dataset}}/download_completed.txt"
    output:
        h5ad=f"{DATA_DIR}/sparse_h5/{{dataset}}_sparse.h5ad"
    shell:
        """
        python {input.script} --dataset {wildcards.dataset} --raw-dir {DATA_DIR}/raw_geo --out-h5 {output.h5ad}
        """

rule preprocess_dataset_to_h5:
    input:
        script="scripts/ingestion/03_preprocess_datasets.py",
        h5ad=f"{DATA_DIR}/sparse_h5/{{dataset}}_sparse.h5ad"
    output:
        h5ad=f"{DATA_DIR}/processed_h5/{{dataset}}_processed.h5ad"
    shell:
        """
        python {input.script} --input-h5 {input.h5ad} --out-h5 {output.h5ad}
        """

rule plot_dataset_embedding:
    input:
        script="scripts/ingestion/04_plot_dataset_embeddings.py",
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
# Sade-Feldman Multi-Resolution Deconvolution, Reference Selection & Validation
# ==============================================================================
SF_REF_DIR = f"{MANUSCRIPT_OUT_DIR}/3f1_build_sf_reference"
INTEG_REF_DIR = f"{MANUSCRIPT_OUT_DIR}/3f2_build_integrated_reference"
DS_REF_DIR = f"{MANUSCRIPT_OUT_DIR}/3f3_build_dataset_reference"
COMB_REF_DIR = f"{MANUSCRIPT_OUT_DIR}/3f4_build_combined_references"
MULTIRES_DECONV_DIR = f"{MANUSCRIPT_OUT_DIR}/3f5_deconvolute_iatlas_multires"
MULTIRES_LR_DIR = f"{MANUSCRIPT_OUT_DIR}/3f6_logistic_regression_multires"
MULTIRES_MILO_DIR = f"{MANUSCRIPT_OUT_DIR}/3f7_analyze_milopy_multires"
MULTIRES_CONC_DIR = f"{MANUSCRIPT_OUT_DIR}/3f8_compare_concordance_multires"
MULTIRES_FIG_DIR = f"{MANUSCRIPT_OUT_DIR}/4d1_plot_sf_multires_validation_figures"
REF_SEL_FIG_DIR = f"{MANUSCRIPT_OUT_DIR}/4d2_plot_cohort_auc_reference_selection"
AUX_REF_DATASETS = ["jerby", "maynard", "ma", "yost"]

rule build_sade_feldman_reference:
    input:
        script="scripts/manuscript_pipeline/3f1_build_sf_reference.py",
        adata=f"{DATA_DIR}/preprocessed/GSE120575.h5ad",
    output:
        phi=f"{SF_REF_DIR}/reference_phi.parquet",
        markers=f"{SF_REF_DIR}/reference_marker_genes.parquet",
        phi_res05=f"{SF_REF_DIR}/reference_phi_res0.5.parquet",
        phi_res10=f"{SF_REF_DIR}/reference_phi_res1.0.parquet",
        phi_res15=f"{SF_REF_DIR}/reference_phi_res1.5.parquet",
        phi_res20=f"{SF_REF_DIR}/reference_phi_res2.0.parquet",
        metrics=f"{SF_REF_DIR}/reference_resolution_metrics_sf.parquet",
    params:
        out_dir=SF_REF_DIR,
    shell:
        """
        python {input.script} --adata {input.adata} --out-dir {params.out_dir}
        """

rule build_integrated_reference:
    input:
        script="scripts/manuscript_pipeline/3f2_build_integrated_reference.py",
        adata=f"{DATA_DIR}/preprocessed/GSE120575.h5ad",
    output:
        phi=f"{INTEG_REF_DIR}/integrated_reference_phi.parquet",
        markers=f"{INTEG_REF_DIR}/integrated_marker_genes.parquet",
        phi_res05=f"{INTEG_REF_DIR}/integrated_reference_phi_res0.5.parquet",
        phi_res10=f"{INTEG_REF_DIR}/integrated_reference_phi_res1.0.parquet",
        phi_res15=f"{INTEG_REF_DIR}/integrated_reference_phi_res1.5.parquet",
        phi_res20=f"{INTEG_REF_DIR}/integrated_reference_phi_res2.0.parquet",
        meta=f"{INTEG_REF_DIR}/integrated_cell_metadata.parquet",
        metrics=f"{INTEG_REF_DIR}/integration_quality_metrics.parquet",
        res_metrics=f"{INTEG_REF_DIR}/reference_resolution_metrics_comb.parquet",
    params:
        out_dir=INTEG_REF_DIR,
    shell:
        """
        python {input.script} --adata-sf {input.adata} --out-dir {params.out_dir}
        """

rule build_dataset_references:
    input:
        script="scripts/manuscript_pipeline/3f3_build_dataset_reference.py",
    output:
        phi_res05=f"{DS_REF_DIR}/reference_phi_{{ref_ds}}_res0.5.parquet",
        metrics=f"{DS_REF_DIR}/reference_resolution_metrics_{{ref_ds}}.parquet",
    wildcard_constraints:
        ref_ds="jerby|maynard|ma|yost"
    params:
        raw_dir=f"{DATA_DIR}/raw",
        out_dir=DS_REF_DIR,
    shell:
        """
        python {input.script} --dataset {wildcards.ref_ds} --raw-dir {params.raw_dir} --out-dir {params.out_dir}
        """

rule build_combined_references:
    input:
        script="scripts/manuscript_pipeline/3f4_build_combined_references.py",
        ds_metrics=expand(f"{DS_REF_DIR}/reference_resolution_metrics_{{ref_ds}}.parquet", ref_ds=AUX_REF_DATASETS),
    output:
        all4_phi=f"{COMB_REF_DIR}/reference_phi_all_4_combined_res0.5.parquet",
        all4_metrics=f"{COMB_REF_DIR}/reference_resolution_metrics_all_4_combined.parquet",
    params:
        raw_dir=f"{DATA_DIR}/raw",
        out_dir=COMB_REF_DIR,
    shell:
        """
        python {input.script} --raw-dir {params.raw_dir} --out-dir {params.out_dir}
        """

rule deconvolute_iatlas_sade_feldman:
    input:
        script="scripts/manuscript_pipeline/3f5_deconvolute_iatlas_multires.py",
        phi=f"{SF_REF_DIR}/reference_phi.parquet",
        integ_phi=f"{INTEG_REF_DIR}/integrated_reference_phi.parquet",
        ds_metrics=expand(f"{DS_REF_DIR}/reference_resolution_metrics_{{ref_ds}}.parquet", ref_ds=AUX_REF_DATASETS),
        all4_phi=f"{COMB_REF_DIR}/reference_phi_all_4_combined_res0.5.parquet",
    output:
        fracs=f"{MULTIRES_DECONV_DIR}/deconv_fractions.parquet",
        fracs_sf_05=f"{MULTIRES_DECONV_DIR}/deconv_fractions_sf_res0.5.parquet",
        fracs_sf_10=f"{MULTIRES_DECONV_DIR}/deconv_fractions_sf_res1.0.parquet",
        fracs_sf_15=f"{MULTIRES_DECONV_DIR}/deconv_fractions_sf_res1.5.parquet",
        fracs_sf_20=f"{MULTIRES_DECONV_DIR}/deconv_fractions_sf_res2.0.parquet",
        fracs_comb_05=f"{MULTIRES_DECONV_DIR}/deconv_fractions_comb_res0.5.parquet",
        fracs_comb_10=f"{MULTIRES_DECONV_DIR}/deconv_fractions_comb_res1.0.parquet",
        fracs_comb_15=f"{MULTIRES_DECONV_DIR}/deconv_fractions_comb_res1.5.parquet",
        fracs_comb_20=f"{MULTIRES_DECONV_DIR}/deconv_fractions_comb_res2.0.parquet",
    params:
        out_dir=MULTIRES_DECONV_DIR,
    shell:
        """
        python {input.script} --out-dir {params.out_dir} --benchmark-mode
        """

rule logistic_regression_sade_feldman:
    input:
        script="scripts/manuscript_pipeline/3f6_logistic_regression_multires.py",
        fracs=f"{MULTIRES_DECONV_DIR}/deconv_fractions.parquet",
        fracs_comb=f"{MULTIRES_DECONV_DIR}/deconv_fractions_comb_res0.5.parquet",
    output:
        results=f"{MULTIRES_LR_DIR}/logistic_regression_results.parquet",
        all_results=f"{MULTIRES_LR_DIR}/logistic_regression_results_all_resolutions.parquet",
    params:
        out_dir=MULTIRES_LR_DIR,
    shell:
        """
        python {input.script} --out-dir {params.out_dir} --benchmark-mode
        """

rule analyze_milopy_sade_feldman:
    input:
        script="scripts/manuscript_pipeline/3f7_analyze_milopy_multires.py",
        adata=f"{DATA_DIR}/preprocessed/GSE120575.h5ad",
    output:
        states=f"{MULTIRES_MILO_DIR}/milopy_cell_state_da.parquet",
        all_states=f"{MULTIRES_MILO_DIR}/milopy_cell_state_da_all.parquet",
        all_res_states=f"{MULTIRES_MILO_DIR}/milopy_cell_state_da_all_resolutions.parquet",
        pre_states=f"{MULTIRES_MILO_DIR}/milopy_cell_state_da_Pre.parquet",
        post_states=f"{MULTIRES_MILO_DIR}/milopy_cell_state_da_Post.parquet",
        cells=f"{MULTIRES_MILO_DIR}/milopy_cell_level_scores.parquet",
    params:
        out_dir=MULTIRES_MILO_DIR,
        milo_dir=f"{DATA_DIR}/output/whole_dataset",
    shell:
        """
        python {input.script} --adata {input.adata} --milo-dir {params.milo_dir} --out-dir {params.out_dir}
        """

rule concordance_sade_feldman:
    input:
        script="scripts/manuscript_pipeline/3f8_compare_concordance_multires.py",
        logistic=f"{MULTIRES_LR_DIR}/logistic_regression_results.parquet",
        logistic_all=f"{MULTIRES_LR_DIR}/logistic_regression_results_all_resolutions.parquet",
        milo_all=f"{MULTIRES_MILO_DIR}/milopy_cell_state_da_all.parquet",
        sf_meta=f"{SF_REF_DIR}/reference_resolution_metrics_sf.parquet",
        comb_meta=f"{INTEG_REF_DIR}/reference_resolution_metrics_comb.parquet",
        ds_metrics=expand(f"{DS_REF_DIR}/reference_resolution_metrics_{{ref_ds}}.parquet", ref_ds=AUX_REF_DATASETS),
        all4_metrics=f"{COMB_REF_DIR}/reference_resolution_metrics_all_4_combined.parquet",
    output:
        metrics=f"{MULTIRES_CONC_DIR}/concordance_metrics.parquet",
        summary=f"{MULTIRES_CONC_DIR}/concordance_summary.parquet",
        cohorts=f"{MULTIRES_CONC_DIR}/cohort_level_concordance_summary.parquet",
        timepoints=f"{MULTIRES_CONC_DIR}/timepoint_concordance_summary.parquet",
        full_metrics=f"{MULTIRES_CONC_DIR}/concordance_metrics_full.parquet",
        benchmark_summary=f"{MULTIRES_CONC_DIR}/multi_resolution_benchmark_summary.parquet",
    params:
        milo_dir=MULTIRES_MILO_DIR,
        out_dir=MULTIRES_CONC_DIR,
    shell:
        """
        python {input.script} --logistic-results {input.logistic} --milo-dir {params.milo_dir} --out-dir {params.out_dir}
        """

rule plot_sade_feldman_validation:
    input:
        script="scripts/manuscript_pipeline/4d1_plot_sf_multires_validation_figures.py",
        markers=f"{SF_REF_DIR}/reference_marker_genes.parquet",
        fracs=f"{MULTIRES_DECONV_DIR}/deconv_fractions.parquet",
        logistic=f"{MULTIRES_LR_DIR}/logistic_regression_results.parquet",
        milo_all=f"{MULTIRES_MILO_DIR}/milopy_cell_state_da_all.parquet",
        cohorts=f"{MULTIRES_CONC_DIR}/cohort_level_concordance_summary.parquet",
        metrics_full=f"{MULTIRES_CONC_DIR}/concordance_metrics_full.parquet",
        cells=f"{MULTIRES_MILO_DIR}/milopy_cell_level_scores.parquet",
        integ_meta=f"{INTEG_REF_DIR}/integrated_cell_metadata.parquet",
        benchmark_summary=f"{MULTIRES_CONC_DIR}/multi_resolution_benchmark_summary.parquet",
    output:
        fig_pred_svg=f"{MULTIRES_FIG_DIR}/fig1_predictive_capacity_benchmark.svg",
        fig_pred_png=f"{MULTIRES_FIG_DIR}/fig1_predictive_capacity_benchmark.png",
        fig_coll_svg=f"{MULTIRES_FIG_DIR}/fig2_signature_collinearity_benchmark.svg",
        fig_gran_svg=f"{MULTIRES_FIG_DIR}/fig3_cellular_granularity_benchmark.svg",
        fig_pareto_svg=f"{MULTIRES_FIG_DIR}/fig4_pareto_tradeoff_optimization.svg",
        fig1=f"{MULTIRES_FIG_DIR}/step01_reference_marker_heatmap.svg",
        fig1b=f"{MULTIRES_FIG_DIR}/step01b_integrated_umap_batch_correction.svg",
        fig1c=f"{MULTIRES_FIG_DIR}/step01c_multi_resolution_marker_heatmaps.svg",
        fig1d=f"{MULTIRES_FIG_DIR}/step01d_combined_reference_umap_resolutions.svg",
        fig2=f"{MULTIRES_FIG_DIR}/step02_deconv_fractions_distribution.svg",
        fig2b=f"{MULTIRES_FIG_DIR}/step02b_cohort_deconv_fractions_grid.svg",
        fig2c=f"{MULTIRES_FIG_DIR}/step02c_cohort_cell_state_stacked_bars.svg",
        fig2d_sf=f"{MULTIRES_FIG_DIR}/step02_deconv_fractions_by_resolution_sf.svg",
        fig2d_comb=f"{MULTIRES_FIG_DIR}/step02_deconv_fractions_by_resolution_comb.svg",
        fig3a=f"{MULTIRES_FIG_DIR}/step03_logistic_regression_volcano.svg",
        fig3b=f"{MULTIRES_FIG_DIR}/step03_logistic_regression_forest.svg",
        fig4a=f"{MULTIRES_FIG_DIR}/step04_milopy_nhood_volcano.svg",
        fig4b=f"{MULTIRES_FIG_DIR}/step04_milopy_cell_state_da.svg",
        fig4c=f"{MULTIRES_FIG_DIR}/step04b_milopy_pre_vs_post.svg",
        fig5=f"{MULTIRES_FIG_DIR}/step05_concordance_scatter.svg",
        fig5b=f"{MULTIRES_FIG_DIR}/step05b_cohort_concordance_comparison.svg",
        fig5c=f"{MULTIRES_FIG_DIR}/step05b_individual_cohort_scatters.svg",
        fig5d_bench=f"{MULTIRES_FIG_DIR}/step05c_resolution_predictive_benchmark.svg",
        fig6=f"{MULTIRES_FIG_DIR}/step06_dual_umap_validation.svg",
        fig6b=f"{MULTIRES_FIG_DIR}/step06b_stratified_umap_validation.svg",
    params:
        data_dir=MULTIRES_CONC_DIR,
        results_dir=MULTIRES_FIG_DIR,
    shell:
        """
        python {input.script} --data-dir {params.data_dir} --results-dir {params.results_dir}
        """

rule plot_cohort_auc_reference_selection:
    input:
        script="scripts/manuscript_pipeline/4d2_plot_cohort_auc_reference_selection.py",
        benchmark_summary=f"{MULTIRES_CONC_DIR}/multi_resolution_benchmark_summary.parquet",
    output:
        svg=f"{REF_SEL_FIG_DIR}/fig1_selected_cohort_auc_and_marginal_gain.svg",
        png=f"{REF_SEL_FIG_DIR}/fig1_selected_cohort_auc_and_marginal_gain.png",
    params:
        data_dir=MULTIRES_CONC_DIR,
        results_dir=REF_SEL_FIG_DIR,
    shell:
        """
        python {input.script} --data-dir {params.data_dir} --results-dir {params.results_dir}
        """

rule run_extended_sade_feldman_pipeline:
    input:
        f"{MULTIRES_FIG_DIR}/fig1_predictive_capacity_benchmark.svg",
        f"{MULTIRES_FIG_DIR}/fig1_predictive_capacity_benchmark.png",
        f"{REF_SEL_FIG_DIR}/fig1_selected_cohort_auc_and_marginal_gain.svg",
        f"{REF_SEL_FIG_DIR}/fig1_selected_cohort_auc_and_marginal_gain.png",
        f"{MULTIRES_FIG_DIR}/step01_reference_marker_heatmap.svg",
        f"{MULTIRES_FIG_DIR}/step01b_integrated_umap_batch_correction.svg",
        f"{MULTIRES_FIG_DIR}/step01c_multi_resolution_marker_heatmaps.svg",
        f"{MULTIRES_FIG_DIR}/step01d_combined_reference_umap_resolutions.svg",
        f"{MULTIRES_FIG_DIR}/step02_deconv_fractions_distribution.svg",
        f"{MULTIRES_FIG_DIR}/step02b_cohort_deconv_fractions_grid.svg",
        f"{MULTIRES_FIG_DIR}/step02c_cohort_cell_state_stacked_bars.svg",
        f"{MULTIRES_FIG_DIR}/step02_deconv_fractions_by_resolution_sf.svg",
        f"{MULTIRES_FIG_DIR}/step02_deconv_fractions_by_resolution_comb.svg",
        f"{MULTIRES_FIG_DIR}/step03_logistic_regression_volcano.svg",
        f"{MULTIRES_FIG_DIR}/step03_logistic_regression_forest.svg",
        f"{MULTIRES_FIG_DIR}/step04_milopy_nhood_volcano.svg",
        f"{MULTIRES_FIG_DIR}/step04_milopy_cell_state_da.svg",
        f"{MULTIRES_FIG_DIR}/step04b_milopy_pre_vs_post.svg",
        f"{MULTIRES_FIG_DIR}/step05_concordance_scatter.svg",
        f"{MULTIRES_FIG_DIR}/step05b_cohort_concordance_comparison.svg",
        f"{MULTIRES_FIG_DIR}/step05b_individual_cohort_scatters.svg",
        f"{MULTIRES_FIG_DIR}/step05c_resolution_predictive_benchmark.svg",
        f"{MULTIRES_FIG_DIR}/step06_dual_umap_validation.svg",
        f"{MULTIRES_FIG_DIR}/step06b_stratified_umap_validation.svg",

rule run_sade_feldman_pipeline:
    input:
        rules.run_extended_sade_feldman_pipeline.input,

# ==============================================================================
# Single-Cell Immuno Datasets Pipeline (Tier 1 Datasets)
# ==============================================================================

rule download_single_cell_immuno_datasets:
    input:
        script="scripts/ingestion/01_download_single_cell_immuno.py"
    output:
        gse120575=f"{DATA_DIR}/raw/GSE120575/GSE120575_tpm.txt.gz",
        gse115978=f"{DATA_DIR}/raw/GSE115978/GSE115978_tpm.csv.gz",
        gse123139=f"{DATA_DIR}/raw/GSE123139/GSE123139_raw_tar.tar",
        gse123813=f"{DATA_DIR}/raw/GSE123813/GSE123813_bcc_counts.txt.gz",
        gse125449=f"{DATA_DIR}/raw/GSE125449/GSE125449_set1_matrix.gz",
        gse179994=f"{DATA_DIR}/raw/GSE179994/GSE179994_counts.gz",
        gse159115=f"{DATA_DIR}/raw/GSE159115/GSE159115_raw_tar.tar",
        gse171306=f"{DATA_DIR}/raw/GSE171306/GSE171306_raw_tar.tar",
        gondal2025=f"{DATA_DIR}/raw/Gondal2025/Gondal2025_zenodo.gz"
    params:
        base_dir=DATA_DIR
    shell:
        """
        set -euo pipefail
        if [ "${{TME_DISABLE_DOWNLOADS:-0}}" = "1" ] || [ -d /ritom/scratch ]; then
            echo "ERROR: Internet dataset downloads are strictly forbidden on CSCS Eiger. Download on olm and transfer via Globus." >&2
            exit 1
        fi
        python {input.script} --base-dir {params.base_dir}
        """


rule preprocess_single_cell_immuno_datasets:
    input:
        script="scripts/ingestion/03_preprocess_datasets.py"
    output:
        marker=f"{DATA_DIR}/preprocessed/preprocess_immuno_completed.txt"
    params:
        input_dir=f"{DATA_DIR}/raw",
        out_dir=f"{DATA_DIR}/preprocessed"
    shell:
        """
        set -euo pipefail
        python {input.script} --input-dir {params.input_dir} --out-dir {params.out_dir}
        touch {output.marker}
        """

rule report_single_cell_immuno_datasets:
    input:
        script="scripts/ingestion/03_generate_single_cell_immuno_report.py"
    output:
        report_md=f"{DATA_DIR}/reports/icb_single_cell_datasets_summary.md",
        svg1=f"{DATA_DIR}/reports/tier1_patient_counts.svg",
        svg2=f"{DATA_DIR}/reports/dataset_tier_summary.svg"
    params:
        base_dir=DATA_DIR
    shell:
        """
        set -euo pipefail
        python {input.script} --base-dir {params.base_dir}
        """

rule check_milopy_compatibility:
    input:
        script="scripts/ingestion/04_check_milopy_compatibility.py",
        marker=f"{DATA_DIR}/preprocessed/preprocess_immuno_completed.txt"
    output:
        report_md=f"{DATA_DIR}/reports/milopy_compatibility_report.md",
        svg_chart=f"{DATA_DIR}/reports/milopy_readiness_summary.svg"
    params:
        base_dir=DATA_DIR
    shell:
        """
        set -euo pipefail
        python {input.script} --base-dir {params.base_dir}
        """

rule run_single_cell_immuno_pipeline:
    input:
        rules.report_single_cell_immuno_datasets.output.report_md,
        rules.report_single_cell_immuno_datasets.output.svg1,
        rules.report_single_cell_immuno_datasets.output.svg2,
        rules.check_milopy_compatibility.output.report_md,
        rules.check_milopy_compatibility.output.svg_chart

rule run_milopy_da:
    input:
        script="scripts/ingestion/05_run_milopy_da.py",
        marker=f"{DATA_DIR}/preprocessed/preprocess_immuno_completed.txt"
    output:
        report_md=f"{DATA_DIR}/reports/milopy_da_cohort_summary.md",
        svg_chart=f"{DATA_DIR}/reports/milopy_da_cohort_summary.svg",
        svg_percentage=f"{DATA_DIR}/reports/milopy_da_cohort_percentage.svg"
    params:
        base_dir=DATA_DIR
    shell:
        """
        set -euo pipefail
        python {input.script} --base-dir {params.base_dir} --dataset all
        """

rule investigate_subtle_datasets:
    input:
        script="scripts/ingestion/06_investigate_subtle_datasets.py",
        marker=f"{DATA_DIR}/preprocessed/preprocess_immuno_completed.txt"
    output:
        report_md=f"{DATA_DIR}/reports/subtle_datasets_investigation.md"
    params:
        base_dir=DATA_DIR
    shell:
        """
        set -euo pipefail
        python {input.script} --base-dir {params.base_dir}
        """

rule integrate_combined_cohorts:
    input:
        script="scripts/ingestion/07_integrate_combined_cohorts.py",
        marker=f"{DATA_DIR}/preprocessed/preprocess_immuno_completed.txt"
    output:
        report_md=f"{DATA_DIR}/reports/integrated_cohorts_summary.md",
        melanoma_da=f"{DATA_DIR}/results/milopy/integrated_melanoma/integrated_da_results.parquet",
        cutaneous_da=f"{DATA_DIR}/results/milopy/integrated_cutaneous/integrated_da_results.parquet",
        pancancer_da=f"{DATA_DIR}/results/milopy/integrated_pancancer/integrated_da_results.parquet"
    params:
        base_dir=DATA_DIR
    shell:
        """
        set -euo pipefail
        python {input.script} --base-dir {params.base_dir} --mode all
        """

MILOPY_COHORTS_9 = [
    "GSE120575",
    "CELLxGENE_7b20c613",
    "CELLxGENE_05a8c945",
    "CELLxGENE_6f9de485",
    "GSE207422",
    "GSE243013",
    "GSE233203",
    "GSE200996",
    "GSE316195",
]
MILOPY_2A1_DIR = f"{MANUSCRIPT_OUT_DIR}/2a1_sample_and_analyze_cohorts"
MILOPY_2A2_DIR = f"{MANUSCRIPT_OUT_DIR}/2a2_combine_cohorts_and_analyze"
MILOPY_4A1_DIR = f"{MANUSCRIPT_OUT_DIR}/4a1_visualize_cohort_umaps"
MILOPY_4A2_DIR = f"{MANUSCRIPT_OUT_DIR}/4a2_generate_combined_figures"

rule sample_and_analyze_milopy_cohorts:
    input:
        script="scripts/manuscript_pipeline/2a1_sample_and_analyze_cohorts.py"
    output:
        da_parquets=expand(f"{MILOPY_2A1_DIR}/{{sc_cohort}}/da_results.parquet", sc_cohort=MILOPY_COHORTS_9),
        cell_parquets=expand(f"{MILOPY_2A1_DIR}/{{sc_cohort}}/cell_level_scores.parquet", sc_cohort=MILOPY_COHORTS_9)
    params:
        base_dir=DATA_DIR,
        out_dir=MILOPY_2A1_DIR,
        cohorts=" ".join(MILOPY_COHORTS_9)
    shell:
        """
        set -euo pipefail
        for cohort in {params.cohorts}; do
            python {input.script} \
                --sample-by custom \
                --value "${{cohort}}" \
                --base-dir {params.base_dir} \
                --out-dir {params.out_dir} \
                --reports-dir {params.out_dir} \
                --no-harmonize
            mkdir -p {params.out_dir}/"${{cohort}}"
            cp -f {params.out_dir}/custom_"${{cohort}}"_none/*.parquet {params.out_dir}/"${{cohort}}"/
            cp -f {params.out_dir}/custom_"${{cohort}}"_none/*.npz {params.out_dir}/"${{cohort}}"/ 2>/dev/null || true
        done
        """

rule visualize_cohort_umaps:
    input:
        script="scripts/manuscript_pipeline/4a1_visualize_cohort_umaps.py",
        da_parquets=expand(f"{MILOPY_2A1_DIR}/{{sc_cohort}}/da_results.parquet", sc_cohort=MILOPY_COHORTS_9),
        cell_parquets=expand(f"{MILOPY_2A1_DIR}/{{sc_cohort}}/cell_level_scores.parquet", sc_cohort=MILOPY_COHORTS_9)
    output:
        completed=f"{MILOPY_4A1_DIR}/visualize_cohort_umaps_completed.txt"
    params:
        base_dir=DATA_DIR,
        results_dir=MILOPY_2A1_DIR,
        reports_dir=MILOPY_4A1_DIR
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --base-dir {params.base_dir} \
            --results-dir {params.results_dir} \
            --reports-dir {params.reports_dir}
        """

rule generate_combined_milopy_figures:
    input:
        script="scripts/manuscript_pipeline/4a2_generate_combined_figures.py",
        da_parquets=expand(f"{MILOPY_2A1_DIR}/{{sc_cohort}}/da_results.parquet", sc_cohort=MILOPY_COHORTS_9),
        cell_parquets=expand(f"{MILOPY_2A1_DIR}/{{sc_cohort}}/cell_level_scores.parquet", sc_cohort=MILOPY_COHORTS_9)
    output:
        volcano_grid=f"{MILOPY_4A2_DIR}/milopy_cross_cohort_volcano_grid.svg",
        umap_grid=f"{MILOPY_4A2_DIR}/milopy_cross_cohort_umaps_logfc.svg",
        volcano_png=f"{MILOPY_4A2_DIR}/milopy_cross_cohort_volcano_grid.png",
        umap_png=f"{MILOPY_4A2_DIR}/milopy_cross_cohort_umaps_logfc.png"
    params:
        results_dir=MILOPY_2A1_DIR,
        reports_dir=MILOPY_4A2_DIR
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --results-dir {params.results_dir} \
            --reports-dir {params.reports_dir} \
            --timepoint all
        """

rule generate_combined_milopy_figures_pre:
    input:
        script="scripts/manuscript_pipeline/4a2_generate_combined_figures.py"
    output:
        volcano_grid=f"{MILOPY_4A2_DIR}/milopy_cross_cohort_volcano_grid_pre.svg",
        umap_grid=f"{MILOPY_4A2_DIR}/milopy_cross_cohort_umaps_logfc_pre.svg",
        volcano_png=f"{MILOPY_4A2_DIR}/milopy_cross_cohort_volcano_grid_pre.png",
        umap_png=f"{MILOPY_4A2_DIR}/milopy_cross_cohort_umaps_logfc_pre.png"
    params:
        results_dir=f"{DATA_DIR}/results/milopy_response_analysis_pre",
        reports_dir=MILOPY_4A2_DIR
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --results-dir {params.results_dir} \
            --reports-dir {params.reports_dir} \
            --timepoint pre
        """

rule generate_combined_milopy_figures_post:
    input:
        script="scripts/manuscript_pipeline/4a2_generate_combined_figures.py"
    output:
        volcano_grid=f"{MILOPY_4A2_DIR}/milopy_cross_cohort_volcano_grid_post.svg",
        umap_grid=f"{MILOPY_4A2_DIR}/milopy_cross_cohort_umaps_logfc_post.svg",
        volcano_png=f"{MILOPY_4A2_DIR}/milopy_cross_cohort_volcano_grid_post.png",
        umap_png=f"{MILOPY_4A2_DIR}/milopy_cross_cohort_umaps_logfc_post.png"
    params:
        results_dir=f"{DATA_DIR}/results/milopy_response_analysis_post",
        reports_dir=MILOPY_4A2_DIR
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --results-dir {params.results_dir} \
            --reports-dir {params.reports_dir} \
            --timepoint post
        """

rule combine_and_run_melanoma_milopy:
    input:
        script="scripts/manuscript_pipeline/2a2_combine_cohorts_and_analyze.py"
    output:
        da_parquet=f"{MILOPY_2A2_DIR}/melanoma/da_results.parquet",
        cell_parquet=f"{MILOPY_2A2_DIR}/melanoma/cell_level_scores.parquet",
        volcano_svg=f"{MILOPY_2A2_DIR}/milopy_combined_melanoma_volcano.svg",
        umap_svg=f"{MILOPY_2A2_DIR}/milopy_combined_melanoma_umap_logfc.svg",
        volcano_png=f"{MILOPY_2A2_DIR}/milopy_combined_melanoma_volcano.png",
        umap_png=f"{MILOPY_2A2_DIR}/milopy_combined_melanoma_umap_logfc.png"
    params:
        group="melanoma",
        base_dir=DATA_DIR,
        out_dir=MILOPY_2A2_DIR,
        reports_dir=MILOPY_2A2_DIR,
        max_cells=25000
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --group {params.group} \
            --base-dir {params.base_dir} \
            --out-dir {params.out_dir} \
            --reports-dir {params.reports_dir} \
            --max-cells-per-cohort {params.max_cells}
        """

rule combine_and_run_nsclc_milopy:
    input:
        script="scripts/manuscript_pipeline/2a2_combine_cohorts_and_analyze.py"
    output:
        da_parquet=f"{MILOPY_2A2_DIR}/nsclc/da_results.parquet",
        cell_parquet=f"{MILOPY_2A2_DIR}/nsclc/cell_level_scores.parquet",
        volcano_svg=f"{MILOPY_2A2_DIR}/milopy_combined_nsclc_volcano.svg",
        umap_svg=f"{MILOPY_2A2_DIR}/milopy_combined_nsclc_umap_logfc.svg",
        volcano_png=f"{MILOPY_2A2_DIR}/milopy_combined_nsclc_volcano.png",
        umap_png=f"{MILOPY_2A2_DIR}/milopy_combined_nsclc_umap_logfc.png"
    params:
        group="nsclc",
        base_dir=DATA_DIR,
        out_dir=MILOPY_2A2_DIR,
        reports_dir=MILOPY_2A2_DIR,
        max_cells=20000
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --group {params.group} \
            --base-dir {params.base_dir} \
            --out-dir {params.out_dir} \
            --reports-dir {params.reports_dir} \
            --max-cells-per-cohort {params.max_cells}
        """

rule combine_and_run_solid_tumors_milopy:
    input:
        script="scripts/manuscript_pipeline/2a2_combine_cohorts_and_analyze.py"
    output:
        da_parquet=f"{MILOPY_2A2_DIR}/solid_tumors/da_results.parquet",
        cell_parquet=f"{MILOPY_2A2_DIR}/solid_tumors/cell_level_scores.parquet",
        volcano_svg=f"{MILOPY_2A2_DIR}/milopy_combined_solid_tumors_volcano.svg",
        umap_svg=f"{MILOPY_2A2_DIR}/milopy_combined_solid_tumors_umap_logfc.svg",
        volcano_png=f"{MILOPY_2A2_DIR}/milopy_combined_solid_tumors_volcano.png",
        umap_png=f"{MILOPY_2A2_DIR}/milopy_combined_solid_tumors_umap_logfc.png"
    params:
        group="solid_tumors",
        base_dir=DATA_DIR,
        out_dir=MILOPY_2A2_DIR,
        reports_dir=MILOPY_2A2_DIR,
        max_cells=20000
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --group {params.group} \
            --base-dir {params.base_dir} \
            --out-dir {params.out_dir} \
            --reports-dir {params.reports_dir} \
            --max-cells-per-cohort {params.max_cells}
        """

rule combine_and_run_pancancer_milopy:
    input:
        script="scripts/manuscript_pipeline/2a2_combine_cohorts_and_analyze.py"
    output:
        da_parquet=f"{MILOPY_2A2_DIR}/pancancer/da_results.parquet",
        cell_parquet=f"{MILOPY_2A2_DIR}/pancancer/cell_level_scores.parquet",
        volcano_svg=f"{MILOPY_2A2_DIR}/milopy_combined_pancancer_volcano.svg",
        umap_svg=f"{MILOPY_2A2_DIR}/milopy_combined_pancancer_umap_logfc.svg",
        volcano_png=f"{MILOPY_2A2_DIR}/milopy_combined_pancancer_volcano.png",
        umap_png=f"{MILOPY_2A2_DIR}/milopy_combined_pancancer_umap_logfc.png"
    params:
        group="pancancer",
        base_dir=DATA_DIR,
        out_dir=MILOPY_2A2_DIR,
        reports_dir=MILOPY_2A2_DIR,
        max_cells=10000
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --group {params.group} \
            --base-dir {params.base_dir} \
            --out-dir {params.out_dir} \
            --reports-dir {params.reports_dir} \
            --max-cells-per-cohort {params.max_cells}
        """

rule combine_and_run_all_milopy:
    input:
        rules.combine_and_run_melanoma_milopy.output.volcano_svg,
        rules.combine_and_run_nsclc_milopy.output.volcano_svg,
        rules.combine_and_run_solid_tumors_milopy.output.volcano_svg,
        rules.combine_and_run_pancancer_milopy.output.volcano_svg

rule visualize_synthetic_umaps:
    input:
        script="scripts/milopy_analysis/visualize_synthetic_umaps.py"
    output:
        parquet="output/synthetic_benchmark/synthetic_cells_umap.parquet",
        composite_svg="results/synthetic_benchmark/synthetic_milo_umaps_composite.svg",
        clusters_svg="results/synthetic_benchmark/synthetic_umap_clusters.svg",
        response_svg="results/synthetic_benchmark/synthetic_umap_response.svg",
        logfc_svg="results/synthetic_benchmark/synthetic_umap_milopy_logfc.svg",
        sig_svg="results/synthetic_benchmark/synthetic_umap_significance.svg"
    params:
        n_cells=3000,
        n_genes=300,
        n_clusters=6,
        n_patients=30,
        fdr=0.10,
        out_dir="output/synthetic_benchmark",
        results_dir="results/synthetic_benchmark"
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --n-cells {params.n_cells} \
            --n-genes {params.n_genes} \
            --n-clusters {params.n_clusters} \
            --n-patients {params.n_patients} \
            --fdr-threshold {params.fdr} \
            --out-dir {params.out_dir} \
            --results-dir {params.results_dir}
        """

rule generate_dataset_tabular_figure:
    input:
        script="scripts/figures/generate_dataset_tabular_figure.py"
    output:
        svg_article="article/figures/dataset_overview/figure_datasets_tabular_compendium.svg",
        png_article="article/figures/dataset_overview/figure_datasets_tabular_compendium.png",
        svg_reports="output/reports/figure_datasets_tabular_compendium.svg",
        png_reports="output/reports/figure_datasets_tabular_compendium.png"
    shell:
        """
        set -euo pipefail
        python {input.script}
        """

rule generate_milopy_method_figure:
    input:
        script="scripts/figures/generate_milopy_workflow.py"
    output:
        svg_article="article/figures/deconvolution/how_milopy_works.svg",
        png_article="article/figures/deconvolution/how_milopy_works.png",
        svg_reports="output/reports/how_milopy_works.svg",
        png_reports="output/reports/how_milopy_works.png"
    shell:
        """
        set -euo pipefail
        python {input.script}
        """

# Single-Cell & Bulk Deconvolution Concordance & Calibration Modules
include: "workflow/rules/concordance.smk"
include: "workflow/rules/calibration.smk"

# Multi-Cohort Stability Selection & Publication Figure Module
include: "workflow/rules/stability_selection.smk"

# Single-Cell Reference Sampling and Classifier HPO Module
include: "workflow/rules/reference_sampling_hpo.smk"

# Complete 37-Script Manuscript Pipeline Module (2a1-4e)
include: "workflow/rules/manuscript_pipeline.smk"

# Consensus-Meta42 vs. Coulton_TAM_Atlas Comparison & ICI Outcome Module
include: "workflow/rules/coulton_tam_comparison.smk"

# Consensus-Meta42 vs. Pan-Cancer T-Cell, CAF & Myeloid Triad Atlases Module
include: "workflow/rules/pan_cancer_atlas_comparison.smk"

# 9-Cohort Single-Cell Reference Deconvolution & Milo Validation Module
include: "workflow/rules/sc9_deconv_validation.smk"








