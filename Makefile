REMOTE_HOST = olm
REMOTE_DIR = ~/python-venv/tme_analysis
SCRIPT = scripts/iAtlas-TMB.py
REMOTE_UV = /home/halu/.local/bin/uv

SCRIPTS := $(wildcard scripts/*.py)
SENTINELS := $(patsubst scripts/%.py,.make-sentinels/%.done,$(SCRIPTS))

.PHONY: sync run-remote run-all-remote

sync:
	ssh $(REMOTE_HOST) "mkdir -p $(REMOTE_DIR)"
	rsync -qavz --delete packages scripts pyproject.toml uv.lock $(REMOTE_HOST):$(REMOTE_DIR)/
	ssh $(REMOTE_HOST) "mkdir -p $(REMOTE_DIR)/output"	
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) sync"

run-remote: sync
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run python $(SCRIPT)"
	mkdir -p output
	rsync -qavz $(REMOTE_HOST):$(REMOTE_DIR)/output/ output/

run-all-remote: sync $(SENTINELS)
	@echo "All changed scripts completed."
	mkdir -p output
	rsync -qavz $(REMOTE_HOST):$(REMOTE_DIR)/output/ output/

.make-sentinels/%.done: scripts/%.py
	@mkdir -p .make-sentinels
	ssh $(REMOTE_HOST) "cd $(REMOTE_DIR) && $(REMOTE_UV) run python scripts/$*.py"
	@touch $@
