REMOTE_HOST = olm
REMOTE_DIR = ~/python-venv/tme_analysis
# TODO: Ensure this matches the DATA_DIR in your Snakefile
REMOTE_DATA_DIR = /storage/halu/data
REMOTE_UV = /home/halu/.local/bin/uv

.PHONY: sync run-remote run-rule run-unified run-sade-feldman download-icb-datasets preprocess-icb-datasets report-icb-datasets check-milopy-compatibility run-milopy-da investigate-subtle-datasets integrate-combined-cohorts run-icb-pipeline pull-results run-all unlock article copy-figures compile-article

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

sync:
	@echo "Syncing code to remote server..."
	ssh $(REMOTE_HOST) "mkdir -p $(REMOTE_DIR)"
	rsync -qavz --delete packages scripts pyproject.toml uv.lock Snakefile $(REMOTE_HOST):$(REMOTE_DIR)/
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) sync"

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
	mkdir -p output/results output/output output/reports
	rsync -qavz $(REMOTE_HOST):$(REMOTE_DATA_DIR)/results/ output/results/
	rsync -qavz $(REMOTE_HOST):$(REMOTE_DATA_DIR)/output/ output/output/
	rsync -qavz $(REMOTE_HOST):$(REMOTE_DATA_DIR)/reports/ output/reports/
