REMOTE_HOST = olm
REMOTE_DIR = ~/python-venv/tme_analysis
# TODO: Ensure this matches the DATA_DIR in your Snakefile
REMOTE_DATA_DIR = /storage/halu/data
REMOTE_UV = /home/halu/.local/bin/uv

.PHONY: sync run-remote pull-results run-all

run-all: sync run-remote pull-results

sync:
	@echo "Syncing code to remote server..."
	ssh $(REMOTE_HOST) "mkdir -p $(REMOTE_DIR)"
	rsync -qavz --delete packages scripts pyproject.toml uv.lock Snakefile $(REMOTE_HOST):$(REMOTE_DIR)/
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) sync"

run-remote:
	@echo "Running Snakemake on remote server..."
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run snakemake --cores all"

pull-results:
	@echo "Pulling results back to local machine..."
	mkdir -p output/results
	rsync -qavz $(REMOTE_HOST):$(REMOTE_DATA_DIR)/results/ output/results/
