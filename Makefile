REMOTE_USER ?= lhaeuser
REMOTE_HOST ?= eiger
SLURM_ACCOUNT ?= em25
SLURM_PARTITION ?= normal
SLURM_STAGE_PARTITION ?= prepost

# Persistent storage ($STORE) - NEVER accessed directly by compute jobs
REMOTE_STORE_BASE = /capstor/store/cscs/empa/em25/halu
REMOTE_STORE_CODE_DIR = $(REMOTE_STORE_BASE)/code
REMOTE_STORE_DATA_DIR = $(REMOTE_STORE_BASE)/data-test

# Active runtime storage ($SCRATCH) - EXCLUSIVE target for all runtime reads and writes
REMOTE_SCRATCH_BASE = /ritom/scratch/cscs/$(REMOTE_USER)
REMOTE_DIR = $(REMOTE_SCRATCH_BASE)/code
REMOTE_DATA_DIR = $(REMOTE_SCRATCH_BASE)/data-test
REMOTE_TEST_DATA_DIR = $(REMOTE_DATA_DIR)
REMOTE_TMP_DIR = $(REMOTE_SCRATCH_BASE)/tmp
REMOTE_UV = /users/$(REMOTE_USER)/.local/bin/uv
LOCAL_UV ?= uv

STAGE_SUBDIR ?=

REMOTE_ENV = TME_DATA_DIR=$(REMOTE_DATA_DIR) \
             TME_DISABLE_DOWNLOADS=1 \
             DISCOVERED_DATA_DIR=$(REMOTE_DATA_DIR) \
             TMPDIR=$(REMOTE_TMP_DIR) \
             UV_CACHE_DIR=$(REMOTE_SCRATCH_BASE)/.cache/uv \
             UV_PYTHON_INSTALL_DIR=$(REMOTE_SCRATCH_BASE)/.local/share/uv/python \
             UV_TOOL_DIR=$(REMOTE_SCRATCH_BASE)/.local/share/uv/tools \
             PYENSEMBL_CACHE_DIR=$(REMOTE_DATA_DIR)/ensembl

SRUN_STAGE = srun -A $(SLURM_ACCOUNT) -p $(SLURM_STAGE_PARTITION) -n 1 -c 8 --mem=16G -t 00:30:00
SRUN_COMPUTE = srun -A $(SLURM_ACCOUNT) -p $(SLURM_PARTITION) -n 1 -c $(WORKERS) --mem=$(MAX_RAM_GB)G -t 06:00:00
SNAKEMAKE_SLURM = $(REMOTE_ENV) $(REMOTE_UV) run snakemake --profile workflow/profiles/eiger

.PHONY: push sync sync-code remote-sync-env remote-test remote-dryrun stage-in stage-out dry-run run-remote run-rule run-unified run-sade-feldman download-icb-datasets preprocess-icb-datasets report-icb-datasets check-milopy-compatibility run-milopy-da analyze-response-cohorts-milopy visualize-cohort-umaps combined-milopy-figures visualize-synthetic-umaps investigate-subtle-datasets integrate-combined-cohorts run-icb-pipeline pull-results run-all unlock article copy-figures compile-article jupyter-start jupyter-status jupyter-stop jupyter-logs jupyter-tunnel jupyter-pull jupyter-kernel setup-manual-downloads setup-manual-downloads-local parallel-download-preprocess remote-preprocess-start remote-preprocess-status remote-preprocess-logs remote-preprocess-stop remote-preprocess-attach download-discovered-test preprocess-discovered-test preprocess-discovered-bg preprocess-discovered-status preprocess-discovered-logs preprocess-discovered-attach preprocess-discovered-stop snakemake-discovered-test stability-selection stability-selection-benchmark stability-selection-plots verify-cibersortx verify-cibersortx-local benchmark-cibersortx benchmark-cibersortx-local hpo hpo-remote hpo-benchmark-remote hpo-local hpo-benchmark-local hpo-cohort hpo-cohort-benchmark-remote hpo-cohort-benchmark-local hpo-all-cohorts hpo-all-cohorts-benchmark-remote hpo-all-cohorts-benchmark-local sample-and-analyze-cohorts sample-and-analyze-cohorts-all plot-milopy-cross-cohort-summary

JUPYTER_PORT ?= 8888
WORKERS ?= 8
MAX_RAM_GB ?= 140
TIER ?= tier0_benchmark
DISCOVERED_RULE ?= process_all_tier0_benchmark_cohorts

# ==============================================================================
# Code Sync & $STORE <-> $SCRATCH Staging (Strictly via Slurm)
# ==============================================================================
sync-code: sync

push:
	@echo "Pushing code to persistent \$$STORE ($(REMOTE_STORE_CODE_DIR)) and active \$$SCRATCH ($(REMOTE_DIR))..."
	ssh $(REMOTE_HOST) "mkdir -p $(REMOTE_STORE_CODE_DIR) $(REMOTE_STORE_CODE_DIR)/data/registry $(REMOTE_DIR) $(REMOTE_DIR)/data/registry $(REMOTE_DATA_DIR) $(REMOTE_TMP_DIR) $(REMOTE_SCRATCH_BASE)/.cache/uv $(REMOTE_SCRATCH_BASE)/slurm_logs"
	rsync -qavz --delete packages scripts config workflow tests .agents pyproject.toml uv.lock Snakefile Makefile $(REMOTE_HOST):$(REMOTE_STORE_CODE_DIR)/
	rsync -qavz data/registry/ $(REMOTE_HOST):$(REMOTE_STORE_CODE_DIR)/data/registry/
	rsync -qavz --update --exclude='data/' jupyter $(REMOTE_HOST):$(REMOTE_STORE_CODE_DIR)/
	rsync -qavz --delete packages scripts config workflow tests .agents pyproject.toml uv.lock Snakefile Makefile $(REMOTE_HOST):$(REMOTE_DIR)/
	rsync -qavz data/registry/ $(REMOTE_HOST):$(REMOTE_DIR)/data/registry/
	rsync -qavz --update --exclude='data/' jupyter $(REMOTE_HOST):$(REMOTE_DIR)/
	@if [ -f .env ]; then \
		echo "Syncing .env credentials to $(REMOTE_HOST)..."; \
		rsync -qavz .env $(REMOTE_HOST):$(REMOTE_STORE_CODE_DIR)/.env; \
		rsync -qavz .env $(REMOTE_HOST):$(REMOTE_DIR)/.env; \
		ssh $(REMOTE_HOST) "chmod 600 $(REMOTE_STORE_CODE_DIR)/.env $(REMOTE_DIR)/.env"; \
	fi

remote-sync-env:
	@echo "Syncing Python environment (.venv) on \$$SCRATCH ($(REMOTE_DIR)) via Slurm..."
	ssh $(REMOTE_HOST) "$(SRUN_STAGE) bash -c 'cd $(REMOTE_DIR) && $(REMOTE_ENV) $(REMOTE_UV) sync'"

sync: push remote-sync-env

remote-test:
	@echo "Running unit tests on $(REMOTE_HOST) \$$SCRATCH via Slurm..."
	ssh $(REMOTE_HOST) "$(SRUN_STAGE) bash -c 'cd $(REMOTE_DIR) && $(REMOTE_ENV) $(REMOTE_UV) run pytest tests/test_ssh_gate.py packages/tme_datasets/tests/test_paths_and_caching.py packages/tme_datasets/tests/test_hpo_cli_script.py packages/single_cell_immuno_datasets/tests/test_package.py -v'"

remote-dryrun:
	@echo "Running Snakemake dry-run on $(REMOTE_HOST) \$$SCRATCH..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_ENV) $(REMOTE_UV) run snakemake $(RULE) --dry-run"

stage-in:
	@echo "Staging data from persistent \$$STORE ($(REMOTE_STORE_DATA_DIR)/$(STAGE_SUBDIR)) to \$$SCRATCH ($(REMOTE_DATA_DIR)/$(STAGE_SUBDIR)) via Slurm..."
	ssh $(REMOTE_HOST) "mkdir -p $(REMOTE_DATA_DIR)/$(STAGE_SUBDIR) && $(SRUN_STAGE) bash -c 'rsync -av --update $(REMOTE_STORE_DATA_DIR)/$(STAGE_SUBDIR) $(REMOTE_DATA_DIR)/$(STAGE_SUBDIR) && ln -sf GSE176249.h5ad $(REMOTE_DATA_DIR)/preprocessed/GSE172433.h5ad'"

stage-out:
	@echo "Staging permanent outputs from \$$SCRATCH ($(REMOTE_DATA_DIR)) back to persistent \$$STORE ($(REMOTE_STORE_DATA_DIR)) via Slurm..."
	ssh $(REMOTE_HOST) "$(SRUN_STAGE) bash -c 'for d in preprocessed manual_download ensembl registry output results reports; do if [ -d \"$(REMOTE_DATA_DIR)/\$$d\" ]; then mkdir -p \"$(REMOTE_STORE_DATA_DIR)/\$$d\" && rsync -av --update \"$(REMOTE_DATA_DIR)/\$$d/\" \"$(REMOTE_STORE_DATA_DIR)/\$$d/\"; fi; done'"

dry-run:
	@if [ -z "$(RULE)" ]; then \
		echo "Error: RULE is not defined. Usage: make dry-run RULE=<rule_name>"; \
		exit 1; \
	fi
	@echo "Running Snakemake dry-run for '$(RULE)' on $(REMOTE_HOST)..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_ENV) $(REMOTE_UV) run snakemake $(RULE) --dry-run"

# ==============================================================================
# Blocked Internet Download Targets on CSCS Eiger
# ==============================================================================
download-icb-datasets download-discovered-test parallel-download-preprocess remote-preprocess-start:
	@echo "=========================================================================================="
	@echo "ERROR: Internet dataset downloads are strictly prohibited on CSCS Alps Eiger."
	@echo "Download and preprocess raw datasets on 'olm' first, then transfer them manually"
	@echo "via Globus to $(REMOTE_STORE_DATA_DIR) and run 'make stage-in'."
	@echo "=========================================================================================="
	@exit 1

remote-preprocess-status:
	@ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && bash scripts/cluster_ops/run_remote_preprocess.sh status"

remote-preprocess-logs:
	@ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && bash scripts/cluster_ops/run_remote_preprocess.sh logs"

remote-preprocess-stop:
	@ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && bash scripts/cluster_ops/run_remote_preprocess.sh stop"

remote-preprocess-attach:
	@ssh -t $(REMOTE_HOST) "cd $(REMOTE_DIR) && bash scripts/cluster_ops/run_remote_preprocess.sh attach"

article:
	@$(MAKE) -C article article

copy-figures:
	@$(MAKE) -C article copy-figures

compile-article:
	@$(MAKE) -C article compile-article

unlock:
	@echo "Unlocking remote Snakemake directory on \$$SCRATCH..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_ENV) $(REMOTE_UV) run snakemake --unlock"

run-all: sync stage-in run-remote stage-out pull-results

preprocess-icb-datasets: sync stage-in
	@echo "Running Tier 1 ICB single-cell preprocessing pipeline on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) preprocess_single_cell_immuno_datasets"
	$(MAKE) stage-out

report-icb-datasets: sync stage-in
	@echo "Generating ICB single-cell summary report on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) report_single_cell_immuno_datasets"
	$(MAKE) stage-out
	$(MAKE) pull-results

check-milopy-compatibility: sync stage-in
	@echo "Verifying milopy compatibility for all preprocessed ICB datasets on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) check_milopy_compatibility --ignore-incomplete"
	$(MAKE) stage-out
	$(MAKE) pull-results

run-milopy-da: sync stage-in
	@echo "Running full milopy DA analysis across all ready ICB datasets on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) run_milopy_da"
	$(MAKE) stage-out
	$(MAKE) pull-results

combine-cohorts-milopy: sync stage-in
	@echo "Running cross-cohort combined milopy analysis on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) combine_and_run_all_milopy"
	$(MAKE) stage-out
	$(MAKE) pull-results

visualize-cohort-umaps: sync stage-in
	@echo "Generating single-cell Milopy logFC UMAPs across all clinical response cohorts on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) visualize_cohort_umaps"
	$(MAKE) stage-out
	$(MAKE) pull-results

combined-milopy-figures: sync stage-in
	@echo "Generating consolidated multi-panel volcano and UMAP figures on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) generate_combined_milopy_figures"
	$(MAKE) stage-out
	$(MAKE) pull-results

analyze-response-cohorts-milopy-pre: sync stage-in
	@echo "Running pre-treatment single-cell milopy analysis on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) run_all_milopy_cohorts_pre generate_combined_milopy_figures_pre"
	$(MAKE) stage-out
	$(MAKE) pull-results

analyze-response-cohorts-milopy-post: sync stage-in
	@echo "Running post-treatment single-cell milopy analysis on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) run_all_milopy_cohorts_post generate_combined_milopy_figures_post"
	$(MAKE) stage-out
	$(MAKE) pull-results

combined-milopy-melanoma: sync stage-in
	@echo "Running harmonized Melanoma Milopy analysis on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) combine_and_run_melanoma_milopy"
	$(MAKE) stage-out
	$(MAKE) pull-results

combined-milopy-nsclc: sync stage-in
	@echo "Running harmonized NSCLC Milopy analysis on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) combine_and_run_nsclc_milopy"
	$(MAKE) stage-out
	$(MAKE) pull-results

combined-milopy-pancancer: sync stage-in
	@echo "Running harmonized Pan-Cancer Milopy analysis on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) combine_and_run_pancancer_milopy"
	$(MAKE) stage-out
	$(MAKE) pull-results

SAMPLE_BY ?= all
SAMPLE_VALUE ?=
SAMPLE_CELLS ?= 15000
SAMPLE_PERTURBATION ?= none
SAMPLE_EXTRA_FLAGS ?=

sample-and-analyze-cohorts: sync stage-in
	@echo "Running multi-cohort Milopy analysis (sample-by: $(SAMPLE_BY)) on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "$(SRUN_COMPUTE) bash -c 'cd $(REMOTE_DIR) && PYTHONPATH=scripts/manuscript_pipeline:\$$PYTHONPATH $(REMOTE_ENV) $(REMOTE_UV) run python scripts/manuscript_pipeline/2a1_sample_and_analyze_cohorts.py \
		--sample-by $(SAMPLE_BY) \
		$(if $(SAMPLE_VALUE),--value $(SAMPLE_VALUE)) \
		--cells-per-dataset $(SAMPLE_CELLS) \
		--perturbation $(SAMPLE_PERTURBATION) \
		--base-dir $(REMOTE_DATA_DIR) \
		--out-dir $(REMOTE_DATA_DIR)/results/milopy_sampled_analysis \
		--reports-dir $(REMOTE_DATA_DIR)/reports/milopy_sampled_analysis \
		$(SAMPLE_EXTRA_FLAGS)'"
	$(MAKE) stage-out
	$(MAKE) pull-results

sample-and-analyze-cohorts-all:
	@$(MAKE) sample-and-analyze-cohorts SAMPLE_BY=all

plot-milopy-cross-cohort-summary:
	@echo "Generating publication-grade portrait cross-cohort differential abundance summary..."
	$(LOCAL_UV) run python scripts/milopy_analysis/plot_cross_cohort_comparison.py \
		--output-svg output/reports/milopy_cohort_comparison_portrait.svg

visualize-synthetic-umaps: sync
	@echo "Running synthetic single-cell UMAP visualizer on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) visualize_synthetic_umaps"
	$(MAKE) stage-out
	$(MAKE) pull-results

dataset-tabular-figure:
	@echo "Generating publication-grade dataset tabular compendium figure..."
	$(LOCAL_UV) run python scripts/figures/generate_dataset_tabular_figure.py

milopy-method-figure:
	@echo "Generating publication-grade Milo method illustration figure..."
	$(LOCAL_UV) run python scripts/figures/generate_milopy_workflow.py

stability-selection:
	@echo "Running stability selection pipeline via Snakemake..."
	$(LOCAL_UV) run snakemake all_stability_selection --cores all

stability-selection-benchmark:
	@echo "Running multi-cohort stability selection benchmark via Snakemake..."
	$(LOCAL_UV) run snakemake run_multi_cohort_stability_benchmark --cores all

stability-selection-plots:
	@echo "Generating stability selection figures via Snakemake..."
	$(LOCAL_UV) run snakemake plot_stability_selection_fitter_overlap --cores all

stability-selection-all-figures:
	@echo "Generating all 4 stability selection publication figures..."
	$(LOCAL_UV) run python scripts/stability_selection_iatlas/plotting/run_all_plots.py --all

investigate-subtle-datasets: sync stage-in
	@echo "Running investigation of subtle datasets GSE123813 and GSE159115 on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) investigate_subtle_datasets"
	$(MAKE) stage-out
	$(MAKE) pull-results

integrate-combined-cohorts: sync stage-in
	@echo "Running combined cohort integration with Harmony on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) integrate_combined_cohorts"
	$(MAKE) stage-out
	$(MAKE) pull-results

run-icb-pipeline: sync stage-in
	@echo "Running full ICB single-cell pipeline on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) run_single_cell_immuno_pipeline"
	$(MAKE) stage-out
	$(MAKE) pull-results

run-unified: sync stage-in
	@echo "Running Snakemake target 'plot_unified_embeddings' on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) plot_unified_embeddings"
	$(MAKE) stage-out

run-sade-feldman: sync stage-in
	@echo "Running Sade-Feldman validation pipeline on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) run_sade_feldman_pipeline"
	$(MAKE) stage-out
	$(MAKE) pull-results

preprocess-discovered-test: sync stage-in
	@echo "Running batch preprocessing for discovered cohorts on $(REMOTE_HOST) via Slurm (RAM cap: 150GB)..."
	ssh $(REMOTE_HOST) "$(SRUN_COMPUTE) bash -c 'mkdir -p $(REMOTE_TEST_DATA_DIR)/preprocessed $(REMOTE_TMP_DIR) && cd $(REMOTE_DIR) && $(REMOTE_ENV) $(REMOTE_UV) run python scripts/manuscript_pipeline/1b_batch_preprocess_cohorts.py --data-raw-dir $(REMOTE_TEST_DATA_DIR)/raw --output-dir $(REMOTE_TEST_DATA_DIR)/preprocessed --max-ram-gb 150.0'"
	$(MAKE) stage-out

preprocess-discovered-bg: sync stage-in
	@echo "Submitting batch preprocessing for discovered cohorts as Slurm batch job on $(REMOTE_HOST)..."
	ssh $(REMOTE_HOST) "mkdir -p $(REMOTE_TEST_DATA_DIR)/preprocessed $(REMOTE_TMP_DIR) $(REMOTE_SCRATCH_BASE)/slurm_logs && \
		sbatch -A $(SLURM_ACCOUNT) -p $(SLURM_PARTITION) -J preprocess_discovered -c $(WORKERS) --mem=150G -t 12:00:00 \
		-o $(REMOTE_SCRATCH_BASE)/slurm_logs/preprocess_discovered.log \
		--wrap='cd $(REMOTE_DIR) && $(REMOTE_ENV) $(REMOTE_UV) run python scripts/manuscript_pipeline/1b_batch_preprocess_cohorts.py --data-raw-dir $(REMOTE_TEST_DATA_DIR)/raw --output-dir $(REMOTE_TEST_DATA_DIR)/preprocessed --max-ram-gb 150.0'"

preprocess-discovered-status:
	@ssh $(REMOTE_HOST) "squeue -u \$$USER -n preprocess_discovered"

preprocess-discovered-logs:
	@ssh $(REMOTE_HOST) "tail -f $(REMOTE_SCRATCH_BASE)/slurm_logs/preprocess_discovered.log"

preprocess-discovered-attach:
	@ssh $(REMOTE_HOST) "squeue -u \$$USER -n preprocess_discovered"

preprocess-discovered-stop:
	@ssh $(REMOTE_HOST) "scancel -u \$$USER -n preprocess_discovered && echo 'Slurm job preprocess_discovered cancelled.'"

snakemake-discovered-test: sync stage-in
	@echo "Running Snakemake workflow for discovered cohorts on $(REMOTE_HOST) via Slurm (target dir: $(REMOTE_TEST_DATA_DIR))..."
	ssh $(REMOTE_HOST) "mkdir -p $(REMOTE_TEST_DATA_DIR)/raw $(REMOTE_TEST_DATA_DIR)/preprocessed $(REMOTE_TMP_DIR) && cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) -s workflow/rules/process_discovered_cohorts.smk $(DISCOVERED_RULE)"
	$(MAKE) stage-out

setup-manual-downloads: sync stage-in
	@echo "Populating manual_download and preprocessed on $(REMOTE_HOST) \$$SCRATCH via Slurm..."
	ssh $(REMOTE_HOST) "$(SRUN_STAGE) bash -c 'cd $(REMOTE_DIR) && $(REMOTE_ENV) $(REMOTE_UV) run python scripts/ingestion/setup_manual_downloads.py'"
	$(MAKE) stage-out

setup-manual-downloads-local:
	@echo "Populating local data/manual_download and data/preprocessed..."
	$(LOCAL_UV) run python scripts/ingestion/setup_manual_downloads.py

jupyter-kernel:
	@echo "Registering tme_analysis ipykernel on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "$(SRUN_STAGE) bash -c 'cd $(REMOTE_DIR) && $(REMOTE_ENV) $(REMOTE_UV) run python -m ipykernel install --user --name tme_analysis --display-name \"Python (tme_analysis)\"'"

jupyter-start: sync jupyter-kernel
	@echo "Starting Jupyter server via Slurm on $(REMOTE_HOST)..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && JUPYTER_PORT=$(JUPYTER_PORT) SLURM_ACCOUNT=$(SLURM_ACCOUNT) bash scripts/cluster_ops/remote_jupyter.sh start"

jupyter-status:
	@ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && JUPYTER_PORT=$(JUPYTER_PORT) bash scripts/cluster_ops/remote_jupyter.sh status"

jupyter-stop:
	@ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && JUPYTER_PORT=$(JUPYTER_PORT) bash scripts/cluster_ops/remote_jupyter.sh stop"

jupyter-logs:
	@ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && bash scripts/cluster_ops/remote_jupyter.sh logs"

jupyter-tunnel:
	@echo "Opening local SSH tunnel to remote Jupyter on port $(JUPYTER_PORT)..."
	@echo "Access at: http://localhost:$(JUPYTER_PORT)/lab"
	@echo "Press Ctrl+C to close the tunnel."
	ssh -N -L $(JUPYTER_PORT):127.0.0.1:$(JUPYTER_PORT) $(REMOTE_HOST)

jupyter-pull:
	@echo "Pulling notebooks from remote server..."
	rsync -qavz --update $(REMOTE_HOST):$(REMOTE_DIR)/jupyter/ jupyter/

run-remote: sync stage-in
	@echo "Running Snakemake on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM)"
	$(MAKE) stage-out

run-rule: sync stage-in
	@if [ -z "$(RULE)" ]; then \
		echo "Error: RULE is not defined. Usage: make run-rule RULE=<rule_name>"; \
		exit 1; \
	fi
	@echo "Running Snakemake rule '$(RULE)' on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) $(RULE)"
	$(MAKE) stage-out

pull-results:
	@echo "Pulling permanent results from $(REMOTE_HOST) \$$STORE ($(REMOTE_STORE_DATA_DIR)) to local machine..."
	mkdir -p output/results output/output output/reports results output
	ssh $(REMOTE_HOST) "mkdir -p $(REMOTE_STORE_DATA_DIR)/results $(REMOTE_STORE_DATA_DIR)/output $(REMOTE_STORE_DATA_DIR)/reports"
	rsync -qavz $(REMOTE_HOST):$(REMOTE_STORE_DATA_DIR)/results/ output/results/ || true
	rsync -qavz $(REMOTE_HOST):$(REMOTE_STORE_DATA_DIR)/output/ output/output/ || true
	rsync -qavz $(REMOTE_HOST):$(REMOTE_STORE_DATA_DIR)/reports/ output/reports/ || true
	rsync -qavz --exclude='.git' $(REMOTE_HOST):$(REMOTE_DIR)/results/ results/ || true
	rsync -qavz --exclude='.git' $(REMOTE_HOST):$(REMOTE_DIR)/output/ output/ || true
	rsync -qavz $(REMOTE_HOST):$(REMOTE_STORE_DATA_DIR)/results/reference_sampling_hpo/ results/reference_sampling_hpo/ || true
	rsync -qavz $(REMOTE_HOST):$(REMOTE_STORE_DATA_DIR)/output/reference_sampling_hpo/ output/reference_sampling_hpo/ || true

verify-cibersortx: sync
	@echo "Running CIBERSORTx verification on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "$(SRUN_STAGE) bash -c 'cd $(REMOTE_DIR) && $(REMOTE_ENV) $(REMOTE_UV) run python scripts/regularized_deconv/verify_cibersortx_docker.py'"

verify-cibersortx-local:
	PYTHONPATH=.venv/lib/python3.14/site-packages:packages/bayesprism/src:packages/instaprism/src \
	/opt/homebrew/bin/python3.14 scripts/regularized_deconv/verify_cibersortx_docker.py

benchmark-cibersortx: sync stage-in
	@echo "Running full CIBERSORTx collinearity benchmark on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "$(SRUN_COMPUTE) bash -c 'cd $(REMOTE_DIR) && $(REMOTE_ENV) $(REMOTE_UV) run python scripts/regularized_deconv/benchmark_collinearity_regularized_deconv.py'"
	$(MAKE) stage-out
	$(MAKE) pull-results

benchmark-cibersortx-local:
	PYTHONPATH=.venv/lib/python3.14/site-packages:packages/bayesprism/src:packages/instaprism/src \
	/opt/homebrew/bin/python3.14 scripts/regularized_deconv/benchmark_collinearity_regularized_deconv.py \
		--cibersortx-username "$(CIBERSORTX_USERNAME)" \
		--cibersortx-token "$(CIBERSORTX_TOKEN)"

plot-deconv-benchmarks:
	@echo "Re-generating all synthetic deconvolution benchmark figures via plotting_utils..."
	PYTHONPATH=.venv/lib/python3.14/site-packages:packages/plotting_utils/src:packages/bayesprism/src:packages/instaprism/src:scripts \
	/opt/homebrew/bin/python3.14 scripts/regularized_deconv/plotting/plot_collinearity_benchmark.py
	PYTHONPATH=.venv/lib/python3.14/site-packages:packages/plotting_utils/src:packages/bayesprism/src:packages/instaprism/src:scripts \
	/opt/homebrew/bin/python3.14 scripts/regularized_deconv/plotting/plot_variance_and_dropout_diagnostics.py
	PYTHONPATH=.venv/lib/python3.14/site-packages:packages/plotting_utils/src:packages/bayesprism/src:packages/instaprism/src:scripts \
	/opt/homebrew/bin/python3.14 scripts/regularized_deconv/plotting/plot_superior_regimes.py
	PYTHONPATH=.venv/lib/python3.14/site-packages:packages/plotting_utils/src:packages/bayesprism/src:packages/instaprism/src:scripts \
	/opt/homebrew/bin/python3.14 scripts/regularized_deconv/plotting/plot_synthetic_overview_and_truth.py
	PYTHONPATH=.venv/lib/python3.14/site-packages:packages/plotting_utils/src:packages/bayesprism/src:packages/instaprism/src:scripts \
	/opt/homebrew/bin/python3.14 scripts/regularized_deconv/plotting/plot_unified_bias_variance_landscape.py
	@echo "All benchmark figures exported cleanly to article/figures/deconvolution/"

benchmark-unified-deconv:
	@echo "Running Unified 3-Factorial Synthetic Deconvolution Benchmark..."
	PYTHONPATH=.venv/lib/python3.14/site-packages:packages/plotting_utils/src:packages/bayesprism/src:packages/instaprism/src:scripts \
	/opt/homebrew/bin/python3.14 -m regularized_deconv.benchmark_unified_synthetic_experiment \
		--collinearities 0.0 0.6 0.9 0.99 \
		--counts 2500 10000 80000 \
		--architectures balanced skewed state_dropout \
		--b-replicates 15 \
		--n-samples 20

# ==============================================================================
# Single-Cell Reference Sampling & Classifier HPO Pipeline (CSCS Slurm Only)
# ==============================================================================
HPO_TRIALS ?= 10
HPO_CANCER ?= melanoma

hpo: hpo-remote

hpo-remote: sync stage-in
	@echo "Running Single-Cell Reference Sampling and Classifier HPO on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) run_reference_sampling_hpo --config hpo_n_trials=$(HPO_TRIALS) hpo_cancer_types=$(HPO_CANCER)"
	$(MAKE) stage-out
	$(MAKE) pull-results

hpo-benchmark-remote: sync
	@echo "Running Single-Cell HPO synthetic benchmark on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) run_reference_sampling_hpo_benchmark"
	$(MAKE) stage-out
	$(MAKE) pull-results

hpo-local:
	@echo "=========================================================================================="
	@echo "ERROR: Local execution of full single-cell reference sampling HPO is strictly prohibited."
	@echo "The pipeline is computationally intensive and requires Slurm execution on '$(REMOTE_HOST)'."
	@echo "To execute on $(REMOTE_HOST):"
	@echo "    make hpo          (or make hpo-remote)"
	@echo "For local unit testing with synthetic cohorts only, run:"
	@echo "    make hpo-benchmark-local"
	@echo "=========================================================================================="
	@exit 1

hpo-benchmark-local:
	@echo "Running lightweight synthetic HPO benchmark locally..."
	$(LOCAL_UV) run python scripts/manuscript_pipeline/3e2_run_hpo_pipeline.py \
		--benchmark-mode \
		--n-trials 3 \
		--out-dir output/manuscript_pipeline/3e2_run_hpo_pipeline/benchmark \
		--results-dir output/manuscript_pipeline/3e2_run_hpo_pipeline/benchmark

# Cohort-Individual HPO (e.g. Gide, 32 cores via Slurm on eiger)
HPO_COHORT ?= Gide

hpo-cohort: sync stage-in
	@echo "Running Cohort-Individual HPO for '$(HPO_COHORT)' on $(REMOTE_HOST) via Slurm (32 cores)..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) run_single_cohort_hpo --config hpo_target_cohort=$(HPO_COHORT) hpo_n_trials=$(HPO_TRIALS) hpo_cancer_types=$(HPO_CANCER) --forcerun run_single_cohort_hpo"
	$(MAKE) stage-out
	$(MAKE) pull-results

hpo-cohort-benchmark-remote: sync
	@echo "Running Cohort-Individual HPO synthetic benchmark for '$(HPO_COHORT)' on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(SNAKEMAKE_SLURM) run_single_cohort_hpo_benchmark --config hpo_target_cohort=$(HPO_COHORT)"
	$(MAKE) stage-out
	$(MAKE) pull-results

hpo-cohort-benchmark-local:
	@echo "Running lightweight synthetic Cohort-Individual HPO benchmark locally for '$(HPO_COHORT)'..."
	$(LOCAL_UV) run python scripts/reference_sampling_hpo/run_single_cohort_hpo.py \
		--target-cohort $(HPO_COHORT) \
		--benchmark-mode

# Cohort-Individual HPO across All Primary Cohorts & Official Benchmark Integration
HPO_ALL_TRIALS ?= 12

hpo-all-cohorts: sync stage-in
	@echo "Running Multi-Cohort Individual HPO & Benchmarking on $(REMOTE_HOST) via Slurm (32 cores)..."
	ssh $(REMOTE_HOST) "srun -A $(SLURM_ACCOUNT) -p $(SLURM_PARTITION) -n 1 -c 32 --mem=64G -t 06:00:00 bash -c 'cd $(REMOTE_DIR) && $(REMOTE_ENV) $(REMOTE_UV) run python scripts/reference_sampling_hpo/run_all_cohorts_hpo_and_benchmarking.py --n-trials $(HPO_ALL_TRIALS) --n-jobs 32'"
	$(MAKE) stage-out
	$(MAKE) pull-results

hpo-all-cohorts-benchmark-remote: sync
	@echo "Running Multi-Cohort Individual HPO benchmark on $(REMOTE_HOST) via Slurm..."
	ssh $(REMOTE_HOST) "$(SRUN_STAGE) bash -c 'cd $(REMOTE_DIR) && $(REMOTE_ENV) $(REMOTE_UV) run python scripts/reference_sampling_hpo/run_all_cohorts_hpo_and_benchmarking.py --benchmark-mode --n-trials 2 --n-jobs 4'"
	$(MAKE) stage-out
	$(MAKE) pull-results

hpo-all-cohorts-benchmark-local:
	@echo "Running lightweight Multi-Cohort Individual HPO benchmark locally..."
	$(LOCAL_UV) run python scripts/reference_sampling_hpo/run_all_cohorts_hpo_and_benchmarking.py \
		--benchmark-mode \
		--force-local \
		--n-trials 2 \
		--n-jobs 4 \
		--cohorts Gide,Hugo \
		--out-dir output/reference_sampling_hpo/benchmark_all \
		--benchmark-dir output/benchmarks/iatlas_roc_auc

# ------------------------------------------------------------------------------
# Single-Cell Pseudobulk Deconvolution vs. Milopy DA Benchmarks
# ------------------------------------------------------------------------------
BENCHMARK_COHORT ?= GSE120575
BENCHMARK_SWEEPS ?= 5

benchmark-pseudobulk-deconv:
	@echo "Running pseudobulk deconvolution benchmark on cohort '$(BENCHMARK_COHORT)'..."
	$(LOCAL_UV) run python scripts/manuscript_pipeline/3a_sade_feldman_pseudobulk_deconv_benchmark.py --cohort $(BENCHMARK_COHORT)

benchmark-pseudobulk-deconv-all:
	@echo "Running pseudobulk deconvolution benchmark across all 9 cohorts..."
	$(LOCAL_UV) run python scripts/manuscript_pipeline/3a_sade_feldman_pseudobulk_deconv_benchmark.py --cohort all

benchmark-pseudobulk-deconv-sweeps:
	@echo "Running $(BENCHMARK_SWEEPS) baseline sweeps of pseudobulk deconvolution benchmark across all 9 cohorts..."
	$(LOCAL_UV) run python scripts/manuscript_pipeline/3a_sade_feldman_pseudobulk_deconv_benchmark.py --cohort all --n-sweeps $(BENCHMARK_SWEEPS)

benchmark-pseudobulk-deconv-remote: sync stage-in
	@echo "Running $(BENCHMARK_SWEEPS)-sweep empirical pseudobulk deconvolution benchmark (Leiden res=2.0) on $(REMOTE_HOST) via Slurm..."
	rm -rf output/manuscript_pipeline/3a_sade_feldman_pseudobulk_deconv_benchmark
	ssh $(REMOTE_HOST) "$(SRUN_COMPUTE) bash -c 'cd $(REMOTE_DIR) && rm -rf output/manuscript_pipeline/3a_sade_feldman_pseudobulk_deconv_benchmark && $(REMOTE_ENV) $(REMOTE_UV) run python scripts/manuscript_pipeline/3a_sade_feldman_pseudobulk_deconv_benchmark.py --cohort all --n-sweeps $(BENCHMARK_SWEEPS) --leiden-resolution 2.0 --adata-dir $(REMOTE_DATA_DIR)/preprocessed --milo-dir $(REMOTE_DATA_DIR)/results/milopy_response_analysis'"
	$(MAKE) stage-out
	$(MAKE) pull-results
	$(MAKE) plot-nine-cohorts-scatter-grid

plot-nine-cohorts-scatter-grid:
	@echo "Generating 3x3 composite scatter grid (Panel c, 130mm) across all 9 cohorts..."
	$(LOCAL_UV) run python scripts/manuscript_pipeline/4a5_plot_cross_cohort_deconv_scatter_grid.py
	$(MAKE) plot-nine-cohorts-correlations-bar

plot-nine-cohorts-correlations-bar:
	@echo "Generating correlation bar chart (Panel d, 50mm) and combined Row 2 composite (180mm)..."
	$(LOCAL_UV) run python scripts/manuscript_pipeline/4a6_plot_cross_cohort_deconv_correlations_bar.py

plot-cross-cohort-umaps-grid:
	@echo "Generating 3x3 single-cell UMAP landscape grid (Panel a, 130mm)..."
	$(LOCAL_UV) run python scripts/manuscript_pipeline/4a3_plot_cross_cohort_umaps_grid.py

plot-cross-cohort-da-bars:
	@echo "Generating differential abundance diverging bar plot (Panel b, 50mm)..."
	$(LOCAL_UV) run python scripts/manuscript_pipeline/4a4_plot_cross_cohort_da_bars.py

plot-immunotherapy-2x2-figure:
	@echo "Generating all 4 panels and master 2x2 composite figure (180mm)..."
	$(LOCAL_UV) run python scripts/manuscript_pipeline/4a3_plot_cross_cohort_umaps_grid.py
	$(LOCAL_UV) run python scripts/manuscript_pipeline/4a4_plot_cross_cohort_da_bars.py
	$(LOCAL_UV) run python scripts/manuscript_pipeline/4a5_plot_cross_cohort_deconv_scatter_grid.py
	$(LOCAL_UV) run python scripts/manuscript_pipeline/4a6_plot_cross_cohort_deconv_correlations_bar.py
	$(LOCAL_UV) run python scripts/manuscript_pipeline/4a7_assemble_immunotherapy_2x2_figure.py

test-pseudobulk-deconv:
	@echo "Running complete test suite for 2x2 immunotherapy response figures..."
	$(LOCAL_UV) run pytest tests/test_immunotherapy_2x2_figure.py tests/test_plot_cross_cohort_deconv_scatter_grid.py tests/test_plot_cross_cohort_deconv_correlations_bar.py tests/test_sade_feldman_pseudobulk_deconv.py -v
