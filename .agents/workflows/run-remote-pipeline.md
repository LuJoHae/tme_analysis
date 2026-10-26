---
name: run-remote-pipeline
description: Safely synchronizes code, stages data to $SCRATCH, runs the Snakemake pipeline via Slurm on CSCS Alps Eiger, stages permanent outputs to $STORE, and pulls results back.
---

# Remote Pipeline Execution Workflow (CSCS Alps Eiger)

When the user asks you to "run the pipeline", "execute the code", or calls `/run-remote-pipeline`, follow these steps to ensure execution occurs strictly via Slurm on CSCS Alps Eiger (`eiger`, account `em25`) with all runtime reads/writes confined to `$SCRATCH` (`/ritom/scratch/cscs/lhaeuser/data-test`).

## Step 1: Pre-flight Verification
1. **Check Connection**: Verify SSH connectivity to `eiger` (`ssh -q -o BatchMode=yes -o ConnectTimeout=5 eiger exit`). If this fails, **STOP** immediately and instruct the user to establish the connection to `eiger`.
2. **Verify Paths**:
   - Persistent `$STORE`: `/capstor/store/cscs/empa/em25/halu/data-test` (data) and `/capstor/store/cscs/empa/em25/halu/code` (code).
   - Runtime `$SCRATCH`: `/ritom/scratch/cscs/lhaeuser/data-test` (data) and `/ritom/scratch/cscs/lhaeuser/code` (code).
   - Ensure `Makefile` and `Snakefile` target `/ritom/scratch/cscs/lhaeuser/data-test` for all runtime I/O.

## Step 2: Synchronize, Stage-In, Execute via Slurm, and Stage-Out
Never run compute commands directly in bash on the login node, and never download datasets from the internet on CSCS. Always execute specific rules via the `Makefile`:

1. Run `make sync` to push local changes to `$STORE/halu/code`, mirror to `$SCRATCH/code`, and update `.venv` on `$SCRATCH` via Slurm.
2. Run `make stage-in` to copy required inputs from `/capstor/store/cscs/empa/em25/halu/data-test` to `/ritom/scratch/cscs/lhaeuser/data-test` via Slurm.
3. Run `make run-rule RULE=<rule_name>` to dispatch Snakemake rules as Slurm jobs (`--profile workflow/profiles/eiger`) on `$SCRATCH`, followed by `make stage-out` to persist outputs back to `/capstor/store/cscs/empa/em25/halu/data-test`.
4. Run `make pull-results` to download generated results to the local `output/results/` directory.

## Step 3: Verify and Report
1. Verify that expected files have been populated in the local `output/results/` directory.
2. Present a concise summary of the Slurm execution and pulled files to the user.
