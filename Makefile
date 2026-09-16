REMOTE_HOST = olm
REMOTE_DIR = ~/python-venv/tme_analysis
# TODO: Ensure this matches the DATA_DIR in your Snakefile
REMOTE_DATA_DIR = /storage/halu/data
REMOTE_UV = /home/halu/.local/bin/uv

.PHONY: sync run-remote run-rule run-unified run-sade-feldman pull-results run-all unlock

unlock:
	@echo "Unlocking remote Snakemake directory..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake --unlock"

run-all: sync run-remote pull-results

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
	mkdir -p output/results output/output
	rsync -qavz $(REMOTE_HOST):$(REMOTE_DATA_DIR)/results/ output/results/
	rsync -qavz $(REMOTE_HOST):$(REMOTE_DATA_DIR)/output/ output/output/
