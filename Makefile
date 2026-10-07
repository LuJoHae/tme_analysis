REMOTE_HOST = olm
REMOTE_DIR = ~/python-venv/tme_analysis
# TODO: Ensure this matches the DATA_DIR in your Snakefile
REMOTE_DATA_DIR = /storage/halu/data-test
REMOTE_TEST_DATA_DIR = /storage/halu/data-test
REMOTE_UV = /home/halu/.local/bin/uv
LOCAL_UV ?= uv


.PHONY: sync sync-code run-remote run-rule run-unified run-sade-feldman download-icb-datasets preprocess-icb-datasets report-icb-datasets check-milopy-compatibility run-milopy-da analyze-response-cohorts-milopy visualize-cohort-umaps combined-milopy-figures visualize-synthetic-umaps investigate-subtle-datasets integrate-combined-cohorts run-icb-pipeline pull-results run-all unlock article copy-figures compile-article jupyter-start jupyter-status jupyter-stop jupyter-logs jupyter-tunnel jupyter-pull jupyter-kernel setup-manual-downloads setup-manual-downloads-local parallel-download-preprocess remote-preprocess-start remote-preprocess-status remote-preprocess-logs remote-preprocess-stop remote-preprocess-attach download-discovered-test preprocess-discovered-test snakemake-discovered-test stability-selection stability-selection-benchmark stability-selection-plots verify-cibersortx verify-cibersortx-local benchmark-cibersortx benchmark-cibersortx-local hpo hpo-remote hpo-benchmark-remote hpo-local hpo-benchmark-local

JUPYTER_PORT ?= 8888
WORKERS ?= 8
MAX_RAM_GB ?= 140
TIER ?= tier0_benchmark
DISCOVERED_RULE ?= process_all_tier0_benchmark_cohorts

remote-preprocess-start: sync
	@echo "Starting detached download & preprocessing daemon on remote $(REMOTE_HOST)..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && WORKERS=$(WORKERS) MAX_RAM_GB=$(MAX_RAM_GB) bash scripts/run_remote_preprocess.sh start"

remote-preprocess-status:
	@ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && bash scripts/run_remote_preprocess.sh status"

remote-preprocess-logs:
	@ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && bash scripts/run_remote_preprocess.sh logs"

remote-preprocess-stop:
	@ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && bash scripts/run_remote_preprocess.sh stop"

remote-preprocess-attach:
	@ssh -t $(REMOTE_HOST) "cd $(REMOTE_DIR) && bash scripts/run_remote_preprocess.sh attach"

parallel-download-preprocess: sync
	@echo "Running parallel download & preprocessing on remote $(REMOTE_HOST)..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && TMPDIR=/storage/halu/tmp $(REMOTE_UV) run python scripts/download_and_preprocess_all.py --workers $(WORKERS) --max-ram-gb $(MAX_RAM_GB)"


article:
	@$(MAKE) -C article article

copy-figures:
	@$(MAKE) -C article copy-figures

compile-article:
	@$(MAKE) -C article compile-article

unlock:
	@echo "Unlocking remote Snakemake directory..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake --unlock"

run-all: sync run-remote pull-results

download-icb-datasets: sync
	@echo "Running Tier 1 ICB single-cell download pipeline on remote server..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake download_single_cell_immuno_datasets --cores all --rerun-incomplete --latency-wait 30"

preprocess-icb-datasets: sync
	@echo "Running Tier 1 ICB single-cell preprocessing pipeline on remote server..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake preprocess_single_cell_immuno_datasets --cores all --rerun-incomplete --latency-wait 30"

report-icb-datasets: sync
	@echo "Generating ICB single-cell summary report on remote server..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake report_single_cell_immuno_datasets --cores all --rerun-incomplete --latency-wait 30"
	$(MAKE) pull-results

check-milopy-compatibility: sync
	@echo "Verifying milopy compatibility for all preprocessed ICB datasets on remote server..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake check_milopy_compatibility --cores all --ignore-incomplete --latency-wait 30"
	$(MAKE) pull-results

run-milopy-da: sync
	@echo "Running full milopy DA analysis across all ready ICB datasets on remote server..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake run_milopy_da --cores all --rerun-incomplete --latency-wait 30"
	$(MAKE) pull-results

combine-cohorts-milopy: sync
	@echo "Running cross-cohort combined milopy analysis on remote $(REMOTE_HOST)..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake combine_and_run_all_milopy --cores all --rerun-incomplete --latency-wait 30"
	$(MAKE) pull-results

visualize-cohort-umaps: sync
	@echo "Generating single-cell Milopy logFC UMAPs across all clinical response cohorts on $(REMOTE_HOST)..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake visualize_cohort_umaps --cores all --rerun-incomplete --latency-wait 30"
	$(MAKE) pull-results

combined-milopy-figures: sync
	@echo "Generating consolidated multi-panel volcano and UMAP figures on $(REMOTE_HOST)..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake generate_combined_milopy_figures --cores all --rerun-incomplete --latency-wait 30"
	$(MAKE) pull-results

analyze-response-cohorts-milopy-pre: sync
	@echo "Running pre-treatment single-cell milopy analysis on remote $(REMOTE_HOST)..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake run_all_milopy_cohorts_pre generate_combined_milopy_figures_pre --cores all --rerun-incomplete --latency-wait 30"
	$(MAKE) pull-results

analyze-response-cohorts-milopy-post: sync
	@echo "Running post-treatment single-cell milopy analysis on remote $(REMOTE_HOST)..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake run_all_milopy_cohorts_post generate_combined_milopy_figures_post --cores all --rerun-incomplete --latency-wait 30"
	$(MAKE) pull-results

combined-milopy-melanoma: sync
	@echo "Running harmonized Melanoma Milopy analysis on remote $(REMOTE_HOST)..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake combine_and_run_melanoma_milopy --cores all --rerun-incomplete --latency-wait 30"
	$(MAKE) pull-results

combined-milopy-nsclc: sync
	@echo "Running harmonized NSCLC Milopy analysis on remote $(REMOTE_HOST)..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake combine_and_run_nsclc_milopy --cores all --rerun-incomplete --latency-wait 30"
	$(MAKE) pull-results

combined-milopy-pancancer: sync
	@echo "Running harmonized Pan-Cancer Milopy analysis on remote $(REMOTE_HOST)..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake combine_and_run_pancancer_milopy --cores all --rerun-incomplete --latency-wait 30"
	$(MAKE) pull-results

visualize-synthetic-umaps: sync
	@echo "Running synthetic single-cell UMAP visualizer on remote $(REMOTE_HOST)..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake visualize_synthetic_umaps --cores all --rerun-incomplete --latency-wait 30"
	$(MAKE) pull-results

dataset-tabular-figure:
	@echo "Generating publication-grade dataset tabular compendium figure..."
	python3 scripts/figures/generate_dataset_tabular_figure.py

milopy-method-figure:
	@echo "Generating publication-grade Milo method illustration figure..."
	python3 scripts/figures/generate_milopy_workflow.py

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

investigate-subtle-datasets: sync
	@echo "Running investigation of subtle datasets GSE123813 and GSE159115 on remote server..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake investigate_subtle_datasets --cores all --rerun-incomplete --latency-wait 30"
	$(MAKE) pull-results

integrate-combined-cohorts: sync
	@echo "Running combined cohort integration with Harmony on remote server..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake integrate_combined_cohorts --cores all --rerun-incomplete --latency-wait 30"
	$(MAKE) pull-results

run-icb-pipeline: sync
	@echo "Running full ICB single-cell pipeline on remote server..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake run_single_cell_immuno_pipeline --cores all --rerun-incomplete --latency-wait 30"
	$(MAKE) pull-results

run-unified: sync
	@echo "Running Snakemake target 'plot_unified_embeddings' on remote server..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake plot_unified_embeddings --cores all"

run-sade-feldman: sync
	@echo "Running Sade-Feldman validation pipeline on remote server..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake run_sade_feldman_pipeline --cores all --rerun-incomplete --latency-wait 30"
	$(MAKE) pull-results

sync-code: sync

sync:
	@echo "Syncing code to remote server..."
	ssh $(REMOTE_HOST) "mkdir -p $(REMOTE_DIR) $(REMOTE_DIR)/data/registry"
	rsync -qavz --delete packages scripts config workflow pyproject.toml uv.lock Snakefile $(REMOTE_HOST):$(REMOTE_DIR)/
	rsync -qavz data/registry/ $(REMOTE_HOST):$(REMOTE_DIR)/data/registry/
	rsync -qavz --update --exclude='data/' jupyter $(REMOTE_HOST):$(REMOTE_DIR)/
	@if [ -f .env ]; then \
		echo "Syncing .env credentials to $(REMOTE_HOST)..."; \
		rsync -qavz .env $(REMOTE_HOST):$(REMOTE_DIR)/.env; \
		ssh $(REMOTE_HOST) "chmod 600 $(REMOTE_DIR)/.env"; \
	fi
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) sync"

download-discovered-test: sync
	@echo "Running parallel download for discovered cohorts on remote $(REMOTE_HOST) to $(REMOTE_TEST_DATA_DIR)..."
	ssh $(REMOTE_HOST) "mkdir -p $(REMOTE_TEST_DATA_DIR)/raw /storage/halu/tmp && cd $(REMOTE_DIR) && TMPDIR=/storage/halu/tmp $(REMOTE_UV) run python scripts/download_cohorts.py --data-raw-dir $(REMOTE_TEST_DATA_DIR)/raw --tier $(TIER) --max-workers $(WORKERS)"

preprocess-discovered-test: sync
	@echo "Running batch preprocessing for discovered cohorts on remote $(REMOTE_HOST) (RAM cap: 150GB)..."
	ssh $(REMOTE_HOST) "mkdir -p $(REMOTE_TEST_DATA_DIR)/preprocessed /storage/halu/tmp && cd $(REMOTE_DIR) && TMPDIR=/storage/halu/tmp $(REMOTE_UV) run python scripts/batch_preprocess_cohorts.py --data-raw-dir $(REMOTE_TEST_DATA_DIR)/raw --output-dir $(REMOTE_TEST_DATA_DIR)/preprocessed --max-ram-gb 150.0"

snakemake-discovered-test: sync
	@echo "Running Snakemake workflow for discovered cohorts on remote $(REMOTE_HOST) (target dir: $(REMOTE_TEST_DATA_DIR))..."
	ssh $(REMOTE_HOST) "mkdir -p $(REMOTE_TEST_DATA_DIR)/raw $(REMOTE_TEST_DATA_DIR)/preprocessed /storage/halu/tmp && cd $(REMOTE_DIR) && TMPDIR=/storage/halu/tmp DISCOVERED_DATA_DIR=$(REMOTE_TEST_DATA_DIR) $(REMOTE_UV) run snakemake -s workflow/rules/process_discovered_cohorts.smk $(DISCOVERED_RULE) --cores $(WORKERS) --resources mem_mb=150000 --rerun-incomplete --latency-wait 30"

setup-manual-downloads: sync
	@echo "Populating data/manual_download and data/preprocessed on remote server..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run python scripts/setup_manual_downloads.py"

setup-manual-downloads-local:
	@echo "Populating local data/manual_download and data/preprocessed..."
	uv run python scripts/setup_manual_downloads.py


jupyter-kernel:
	@echo "Registering tme_analysis ipykernel on remote server..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run python -m ipykernel install --user --name tme_analysis --display-name 'Python (tme_analysis)'"

jupyter-start: sync jupyter-kernel
	@echo "Starting Jupyter server on remote server..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && JUPYTER_PORT=$(JUPYTER_PORT) bash scripts/remote_jupyter.sh start"

jupyter-status:
	@ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && JUPYTER_PORT=$(JUPYTER_PORT) bash scripts/remote_jupyter.sh status"

jupyter-stop:
	@ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && JUPYTER_PORT=$(JUPYTER_PORT) bash scripts/remote_jupyter.sh stop"

jupyter-logs:
	@ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && bash scripts/remote_jupyter.sh logs"

jupyter-tunnel:
	@echo "Opening local SSH tunnel to remote Jupyter on port $(JUPYTER_PORT)..."
	@echo "Access at: http://localhost:$(JUPYTER_PORT)/lab"
	@echo "Press Ctrl+C to close the tunnel."
	ssh -N -L $(JUPYTER_PORT):127.0.0.1:$(JUPYTER_PORT) $(REMOTE_HOST)

jupyter-pull:
	@echo "Pulling notebooks from remote server..."
	rsync -qavz --update $(REMOTE_HOST):$(REMOTE_DIR)/jupyter/ jupyter/

run-remote:
	@echo "Running Snakemake on remote server..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake --cores all"

run-rule:
	@if [ -z "$(RULE)" ]; then \
		echo "Error: RULE is not defined. Usage: make run-rule RULE=<rule_name>"; \
		exit 1; \
	fi
	@echo "Running Snakemake rule '$(RULE)' on remote server..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake $(RULE) --cores all --rerun-incomplete --latency-wait 30"

pull-results:
	@echo "Pulling results back to local machine..."
	mkdir -p output/results output/output output/reports results output
	ssh $(REMOTE_HOST) "mkdir -p $(REMOTE_DATA_DIR)/results $(REMOTE_DATA_DIR)/output $(REMOTE_DATA_DIR)/reports"
	rsync -qavz $(REMOTE_HOST):$(REMOTE_DATA_DIR)/results/ output/results/ || true
	rsync -qavz $(REMOTE_HOST):$(REMOTE_DATA_DIR)/output/ output/output/ || true
	rsync -qavz $(REMOTE_HOST):$(REMOTE_DATA_DIR)/reports/ output/reports/ || true
	rsync -qavz --exclude='.git' $(REMOTE_HOST):$(REMOTE_DIR)/results/ results/ || true
	rsync -qavz --exclude='.git' $(REMOTE_HOST):$(REMOTE_DIR)/output/ output/ || true
	rsync -qavz $(REMOTE_HOST):$(REMOTE_DATA_DIR)/results/reference_sampling_hpo/ results/reference_sampling_hpo/ || true
	rsync -qavz $(REMOTE_HOST):$(REMOTE_DATA_DIR)/output/reference_sampling_hpo/ output/reference_sampling_hpo/ || true
verify-cibersortx: sync
	@echo "Running CIBERSORTx Docker verification on remote $(REMOTE_HOST)..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run python scripts/regularized_deconv/verify_cibersortx_docker.py"

verify-cibersortx-local:
	PYTHONPATH=.venv/lib/python3.14/site-packages:packages/bayesprism/src:packages/instaprism/src \
	/opt/homebrew/bin/python3.14 scripts/regularized_deconv/verify_cibersortx_docker.py

benchmark-cibersortx: sync
	@echo "Running full CIBERSORTx collinearity benchmark on remote $(REMOTE_HOST)..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && TMPDIR=/storage/halu/tmp $(REMOTE_UV) run python scripts/regularized_deconv/benchmark_collinearity_regularized_deconv.py"
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
# Single-Cell Reference Sampling & Classifier HPO Pipeline (olm-only)
# ==============================================================================
HPO_TRIALS ?= 10
HPO_CANCER ?= melanoma

hpo: hpo-remote

hpo-remote: sync
	@echo "Running Single-Cell Reference Sampling and Classifier HPO on remote $(REMOTE_HOST)..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake run_reference_sampling_hpo --config hpo_n_trials=$(HPO_TRIALS) hpo_cancer_types=$(HPO_CANCER) --cores all --rerun-incomplete --latency-wait 30"
	$(MAKE) pull-results


hpo-benchmark-remote: sync
	@echo "Running Single-Cell HPO synthetic benchmark on remote $(REMOTE_HOST)..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake run_reference_sampling_hpo_benchmark --cores 2 --rerun-incomplete --latency-wait 30"
	$(MAKE) pull-results

hpo-local:
	@echo "=========================================================================================="
	@echo "ERROR: Local execution of full single-cell reference sampling HPO is strictly prohibited."
	@echo "The pipeline is computationally intensive and requires >64GB RAM on remote host '$(REMOTE_HOST)'."
	@echo "To execute on $(REMOTE_HOST):"
	@echo "    make hpo          (or make hpo-remote)"
	@echo "For local unit testing with synthetic cohorts only, run:"
	@echo "    make hpo-benchmark-local"
	@echo "=========================================================================================="
	@exit 1

hpo-benchmark-local:
	@echo "Running lightweight synthetic HPO benchmark locally..."
	$(LOCAL_UV) run python scripts/reference_sampling_hpo/run_hpo_pipeline.py \
		--benchmark-mode \
		--n-trials 3 \
		--out-dir output/reference_sampling_hpo/benchmark \
		--results-dir results/reference_sampling_hpo/benchmark

